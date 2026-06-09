# -*- coding: utf-8 -*-
# env/geo3D_sim04
#########################

# author: arkriger - 2023 - 2026
# github: https://github.com/AdrianKriger/geo3DopenSim
#########################

import os
import re
import math
import numpy as np
import pandas as pd
from shapely.geometry import Point

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection
from IPython.display import IFrame, display

def load_openfoam_vtk(file_name='hubHeight.vtk', wind_deg=0, extent=None, radius=400.0):
    """
    Generic, state-aware parser for OpenFOAM ASCII VTK files containing CELL_DATA fields.
    Automatically discovers point populations, arbitrary topology shapes, and dynamic FIELD sets.
    """
    vtk_path = file_name #os.path.join(case_path, file_name)
    if not os.path.exists(vtk_path):
        raise FileNotFoundError(f"Could not find VTK file at: {vtk_path}")
        
    #print(f"Executing generic VTK stream analysis on: {file_name}")
    
    # Storage blocks for parsed token arrays
    raw_coords = None
    poly_ints = None
    n_cells = 0
    fields_dict = {}

    # --- Phase 1: Robust State-Aware Linear Scanner ---
    with open(vtk_path, 'r', encoding='utf-8', errors='ignore') as f:
        while True:
            line = f.readline()
            if not line:
                break
                
            line_strip = line.strip()
            if not line_strip:
                continue

            # A. Catch Points Array Layout
            if line_strip.startswith('POINTS'):
                parts = line_strip.split()
                n_points = int(parts[1])
                dtype_str = parts[2].lower()
                #print(f" -> Extracting {n_points} vertices ({dtype_str})...")
                
                # Dynamic stream sizing: 3 coordinates per point
                raw_coords = np.fromfile(f, dtype=float, count=n_points * 3, sep=' ').reshape(n_points, 3)

            # B. Catch Polygon Topology Array Layout
            elif line_strip.startswith('POLYGONS'):
                parts = line_strip.split()
                n_cells = int(parts[1])
                total_ints = int(parts[2])
                #print(f" -> Extracting {n_cells} polygon faces (Descriptor space: {total_ints} ints)...")
                
                poly_ints = np.fromfile(f, dtype=int, count=total_ints, sep=' ')

            # C. Catch Variable Field Block Metadata
            elif line_strip.startswith('CELL_DATA'):
                n_cells_verify = int(line_strip.split()[1])
                assert n_cells == n_cells_verify, "Topology count mismatch encountered during streaming section pass."
                
                # Check next marker line to resolve inner FIELD tags
                next_line = f.readline().strip()
                if next_line.startswith('FIELD'):
                    field_parts = next_line.split()
                    num_fields = int(field_parts[2])
                    
                    # Dynamically ingest every field written inside this block
                    for _ in range(num_fields):
                        attr_meta = f.readline().strip().split()
                        attr_name = attr_meta[0]
                        num_comps = int(attr_meta[1])
                        num_tuples = int(attr_meta[2])
                        # Array data stream type flag
                        data_type = attr_meta[3].lower() 
                        
                        total_elements = num_tuples * num_comps
                        #print(f"   + Discovered dynamic array attribute: '{attr_name}' [{num_comps} component(s), {num_tuples} items]")
                        
                        # Ingest data elements straight into flat numpy structures
                        field_data = np.fromfile(f, dtype=float, count=total_elements, sep=' ')
                        
                        if num_comps > 1:
                            fields_dict[attr_name] = field_data.reshape(num_tuples, num_comps)
                        else:
                            fields_dict[attr_name] = field_data

    # --- Phase 2: Vectorized Spatial Inversion and Polygon Centroid Reconstruction ---
    #print(" -> Processing element cell connectivity midpoints...")
    cell_centers = np.empty((n_cells, 3), dtype=float)
    
    # Rapid pointer navigation tracking across varying layout structures (triangles vs quads)
    idx = 0
    cell_idx = 0
    while idx < len(poly_ints):
        n_verts = poly_ints[idx]
        vert_indices = poly_ints[idx + 1 : idx + 1 + n_verts]
        # Calculate localized geometric midpoint coordinate
        cell_centers[cell_idx] = raw_coords[vert_indices].mean(axis=0)
        
        idx += 1 + n_verts
        cell_idx += 1

    # --- Phase 3: Spatial Transformations & Angle Corrections ---
    # Unroll rotation transformation matrices natively around local mesh origin
    inv_theta = math.radians(270 - wind_deg)
    cos_a, sin_a = np.cos(inv_theta), np.sin(inv_theta)
    
    x_unrot = cell_centers[:, 0] * cos_a - cell_centers[:, 1] * sin_a
    y_unrot = cell_centers[:, 0] * sin_a + cell_centers[:, 1] * cos_a

    # Resolve global anchor coordinate offset shifting
    #if extent is not None:
    x_off = (extent[0] + extent[2]) / 2.0
    y_off = (extent[1] + extent[3]) / 2.0

    final_x = x_unrot + x_off
    final_y = y_unrot + y_off
    final_z = cell_centers[:, 2]

    # --- Phase 4: Dynamic DataFrame Assembly & Masking ---
    dist_sq = (final_x - x_off)**2 + (final_y - y_off)**2
    mask = dist_sq <= radius**2
    
    # Establish base data storage layout
    df_dict = {
        'X': final_x[mask],
        'Y': final_y[mask],
        'Z': final_z[mask]
    }
    
    # Dynamically inject fields parsed from the file metadata
    if 'U' in fields_dict:
        u_raw = fields_dict['U']
        # Realign raw simulation velocity vectors into geographic tracking directions
        df_dict['U'] = (u_raw[:, 0] * cos_a - u_raw[:, 1] * sin_a)[mask]
        df_dict['V'] = (u_raw[:, 0] * sin_a + u_raw[:, 1] * cos_a)[mask]
        df_dict['u_mag'] = np.linalg.norm(u_raw, axis=1)[mask]
    
    if 'k' in fields_dict:
        df_dict['k'] = fields_dict['k'][mask]
        
    # Append any remaining field arrays written into the source (e.g. p, nut, etc.)
    for key, data in fields_dict.items():
        if key not in ['U', 'k'] and len(data) == n_cells:
            df_dict[key] = data[mask]

    df = pd.DataFrame(df_dict)
    #print(f"Completed structural generation. Exported {len(df)} cells inside active bounds.\n")
    
    return df, x_off, y_off

    
def plot_wind_analysis(ped_df, gdf, center_x_utm, center_y_utm, radius=400, title_suffix="xxxx"):
    """
    Generates a side-by-side Velocity Magnitude and Vector Flow plot.
    center_coords: tuple (center_x_utm, center_y_utm)
    """
    #if len(ped_df) > 20000:
    #    bins = 100
    #if len(ped_df) < 10000:
    #    bins = 10
    
    # Create bins
    #ped_df['x_bin'] = (ped_df['X'] // bins) * bins
    #ped_df['y_bin'] = (ped_df['Y'] // bins) * bins

    # Aggregate: Mean velocity per bin
    #binned_df = ped_df.groupby(['x_bin', 'y_bin']).agg({
    #    'U': 'mean', 
    #    'V': 'mean', 
    #    'u_mag': 'mean'
    #}).reset_index()
    
    # 1. Setup the figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 9), sharey=True)

    # --- MAP A: Magnitude (Tricontour) ---
    CNTRsampled_df = ped_df.iloc[::].copy()

    #cntr = ax1.tricontourf(binned_df['x_bin'], binned_df['y_bin'], binned_df['u_mag'],  levels=20, cmap='jet', alpha=0.7)
    cntr = ax1.tricontourf(CNTRsampled_df['X'], CNTRsampled_df['Y'], CNTRsampled_df['u_mag'], cmap='jet', alpha=0.6) #scale=120, alpha=0.9, width=0.003
    ax1.scatter(gdf.geometry.x, gdf.geometry.y, marker='x', color='black', s=20)#, label='Turbines')
    ax1.set_title('A: Wind Velocity Magnitude (m/s)', loc='left', pad=15, weight='bold')

    # --- MAP B: Flow (Quiver) ---
    QUIVsampled_df = ped_df.iloc[::20].copy()
    # Plot binned_df instead of the raw xx rows
    #qv = ax2.quiver(binned_df['x_bin'], binned_df['y_bin'], binned_df['U'], binned_df['V'], binned_df['u_mag'], cmap='jet', scale=120, alpha=0.9, width=0.003)  
    qv = ax2.quiver(QUIVsampled_df['X'], QUIVsampled_df['Y'], QUIVsampled_df['U'], QUIVsampled_df['V'], QUIVsampled_df['u_mag'], cmap='jet', scale=120, alpha=0.9, width=0.003)
    plt.scatter(gdf.geometry.x, gdf.geometry.y, marker='x', color='black', s=20, label='Turbines')

    # --- ADD LABELS TO THE MASTS ON AX2 ---
    for _, row in gdf.iterrows():
        mx = row.geometry.x
        my = row.geometry.y
        m_id = str(row['mastID']).zfill(2) # Formats 1 to '01', 15 to '15'
        
        # Add text slightly offset to the upper-right of the 'x' marker so they don't overlap
        ax2.text(
            mx + 15, my + 15,            # 15-meter coordinate padding shift
            f"B{m_id}",                 # Text string (e.g., 'B15')
            color='black', 
            fontsize=10, 
            weight='bold',
            bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=1) # Clean backdrop
        )

    ax2.set_title('B: Vector Flow Field', loc='left', pad=15, weight='bold')

    # --- AXIS & SPATIAL STYLING ---
    for ax in [ax1, ax2]:
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel('Easting (m)')
        # Force strict 400m AOI
        ax.set_xlim(center_x_utm - radius, center_x_utm + radius)
        ax.set_ylim(center_y_utm - radius, center_y_utm + radius)

    ax1.set_ylabel('Northing (m)')

    # --- SHARED COLORBAR (The "No Squish" Fix) ---
    fig.subplots_adjust(right=0.9) 
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7]) 
    fig.colorbar(cntr, cax=cbar_ax, label='Wind Speed (m/s)')

    plt.suptitle(f'Gouda Wind Facility, Western Cape | SSW @ 600 hour/year - {title_suffix}', fontsize=16, y=0.98)
    
    return fig, (ax1, ax2)
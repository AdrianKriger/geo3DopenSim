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

def reconstruct_openfoam_results(case_path, gdf, wind_deg, extent, radius=400.0):

    def parse_vector(path):
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        match = re.search(r'\n(\d+)\s*\n\s*\(', content)
        n = int(match.group(1))
        limit = content.find('boundaryField', match.end())
        if limit == -1: limit = len(content)
        end_pos = content.rfind(')', match.end(), limit)
        data_str = content[match.end():end_pos]
        return np.fromstring(data_str.replace('(', ' ').replace(')', ' '),
                             dtype=float, sep=' ').reshape(n, 3)

    def parse_labels(path):
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        match = re.search(r'\n(\d+)\s*\n\s*\(', content)
        n = int(match.group(1))
        data_str = content[match.end() : content.rfind(')', match.end())]
        return np.fromstring(data_str, dtype=int, sep=' ')

    def parse_faces_vectorized(path):
        """Returns (face_point_indices, face_sizes) if mixed polyhedral,
           or a single (N_faces, K) array if uniform."""
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        # Extract all face definitions
        raw = re.findall(r'\d+\(([^)]+)\)', content)
        sizes = []
        indices = []
        for face_str in raw:
            pts_idx = np.fromstring(face_str, dtype=np.int32, sep=' ')
            sizes.append(len(pts_idx))
            indices.append(pts_idx)
        return indices, np.array(sizes)

    # Load
    pts     = parse_vector(os.path.join(case_path, 'points'))
    owner   = parse_labels(os.path.join(case_path, 'owner'))
    neighbour = parse_labels(os.path.join(case_path, 'neighbour'))
    u_field = parse_vector(os.path.join(case_path, 'U'))
    faces, face_sizes = parse_faces_vectorized(os.path.join(case_path, 'faces'))

    n_cells = u_field.shape[0]

    # --- Vectorized face centers ---
    # If all faces have the same number of vertices (common in hex meshes), use a single array op
    if np.all(face_sizes == face_sizes[0]):
        k = face_sizes[0]
        flat_idx = np.concatenate(faces).reshape(-1, k)
        face_centers = pts[flat_idx].mean(axis=1)        # (N_faces, 3) — fully vectorized
    else:
        # Mixed polyhedral: still faster than per-face np.mean with Python loop
        flat_idx = np.concatenate(faces)
        offsets   = np.concatenate([[0], np.cumsum(face_sizes)])
        # Use reduceat for a single-pass summation
        sums = np.add.reduceat(pts[flat_idx], offsets[:-1], axis=0)
        face_centers = sums / face_sizes[:, None]

    # --- Fast cell center accumulation with np.bincount ---
    all_cells  = np.concatenate([owner, neighbour])
    all_fcenters = np.vstack([face_centers[np.arange(len(owner))],
                               face_centers[np.arange(len(neighbour))]])
    cell_cx = np.bincount(all_cells, weights=all_fcenters[:, 0], minlength=n_cells)
    cell_cy = np.bincount(all_cells, weights=all_fcenters[:, 1], minlength=n_cells)
    cell_cz = np.bincount(all_cells, weights=all_fcenters[:, 2], minlength=n_cells)
    face_counts = np.bincount(all_cells, minlength=n_cells)

    cell_coords = np.stack([cell_cx, cell_cy, cell_cz], axis=1) / face_counts[:, None]

    x_off = (extent[0] + extent[2]) / 2.0
    y_off = (extent[1] + extent[3]) / 2.0
    
    # Apply inverse rotation to simulation coordinates
    rot_rad =  np.radians(270 - wind_deg)                       
    theta = rot_rad
    cos_i = np.cos(theta)
    sin_i = np.sin(theta)

    x_local_unrot = cell_coords[:, 0] * cos_i - cell_coords[:, 1] * sin_i
    y_local_unrot = cell_coords[:, 0] * sin_i + cell_coords[:, 1] * cos_i

    # 4. Transform back to UTM space by adding the original offsets
    final_x = x_local_unrot + x_off
    final_y = y_local_unrot + y_off
    final_z = cell_coords[:, 2]

    # 5. Rotate the Velocity Vectors back to UTM (North-up)
    # Vectors only need the inverse rotation, no translation
    u_final = u_field[:, 0] * cos_i - u_field[:, 1] * sin_i
    v_final = u_field[:, 0] * sin_i + u_field[:, 1] * cos_i

    # GIS alignment
    cx_utm, cy_utm = gdf.geometry.x.mean(), gdf.geometry.y.mean()

    # Spatial filter
    dist_sq = (final_x - cx_utm)**2 + (final_y - cy_utm)**2
    mask = dist_sq <= radius**2
    
    df = pd.DataFrame({
        'X': final_x[mask], 'Y': final_y[mask], 'Z': final_z[mask],
        #'U': u_field[mask, 0], 
        #'V': u_field[mask, 1],
        'U': u_final[mask], 
        'V': v_final[mask],
        #'u_mag': np.linalg.norm(u_field[mask], axis=1)
        #'u_mag': np.sqrt(u_final[mask]**2 + v_final[mask]**2+ u_field[mask, 2]**2)
        'u_mag': np.sqrt(u_final[mask]**2 + v_final[mask]**2)

    })
    return df, cx_utm, cy_utm

    
def plot_wind_analysis(ped_df, gdf, center_x_utm, center_y_utm, radius=400, title_suffix="xxxx"):
    """
    Generates a side-by-side Velocity Magnitude and Vector Flow plot.
    center_coords: tuple (center_x_utm, center_y_utm)
    """
    if len(ped_df) > 20000:
        bins = 100
    if len(ped_df) < 10000:
        bins = 10
    
    # Create bins
    ped_df['x_bin'] = (ped_df['X'] // bins) * bins
    ped_df['y_bin'] = (ped_df['Y'] // bins) * bins

    # Aggregate: Mean velocity per bin
    binned_df = ped_df.groupby(['x_bin', 'y_bin']).agg({
        'U': 'mean', 
        'V': 'mean', 
        'u_mag': 'mean'
    }).reset_index()
    
    # 1. Setup the figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 9), sharey=True)

    # --- MAP A: Magnitude (Tricontour) ---
    
    cntr = ax1.tricontourf(binned_df['x_bin'], binned_df['y_bin'], binned_df['u_mag'],  levels=20, cmap='jet', alpha=0.7)
    plt.scatter(gdf.geometry.x, gdf.geometry.y, marker='x', color='black', s=20, label='Turbines')
    ax1.set_title('A: Wind Velocity Magnitude (m/s)', loc='left', pad=15, weight='bold')

    # --- MAP B: Flow (Quiver) ---
    
    # Plot binned_df instead of the raw xx rows
    qv = ax2.quiver(binned_df['x_bin'], binned_df['y_bin'], binned_df['U'], binned_df['V'], binned_df['u_mag'], cmap='jet', scale=120, alpha=0.9, width=0.003)    
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
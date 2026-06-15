# -*- coding: utf-8 -*-
# env/geo3DopenSim02
#########################
# author: arkriger - 2023 - 2026
# github: https://github.com/AdrianKriger/geo3D
# Consolidated helper functions for geo3DopenSim. decoupled UTCI - scence.obj mesh with terrain, recreation, roads and buildings
#########################

import os
import re
import datetime
from datetime import timedelta

import math
import numpy as np
import pandas as pd

import shapely 
from shapely.geometry import Point
from shapely.affinity import translate, scale
from shapely.ops import snap, transform, unary_union

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.patches import ConnectionPatch, Rectangle
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
    inv_theta = np.radians(270 - wind_deg)
    cos_a, sin_a = np.cos(inv_theta), np.sin(inv_theta)
    
    x_unrot = cell_centers[:, 0] * cos_a - cell_centers[:, 1] * sin_a
    y_unrot = cell_centers[:, 0] * sin_a + cell_centers[:, 1] * cos_a

    # Resolve global anchor coordinate offset shifting
    if extent is not None:
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
    
    return df#, x_off, y_off


def plot_wind_analysis(ped_df, gdf, center_x_utm, center_y_utm, radius=400, title_suffix="xxxx"):
    """
    Generates a side-by-side Velocity Magnitude and Vector Flow plot.
    center_coords: tuple (center_x_utm, center_y_utm)
    """
    CNTRsampled_df = ped_df.iloc[::10].copy()

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
    
    #cntr = ax1.tricontourf(binned_df['x_bin'], binned_df['y_bin'], binned_df['u_mag'],  levels=20, cmap='jet', alpha=0.7)
    #gdf.plot(ax=ax1, facecolor='grey', edgecolor='grey', alpha=0.1, zorder=1)

    cntr = ax1.tricontourf(CNTRsampled_df['X'], CNTRsampled_df['Y'], CNTRsampled_df['u_mag'], cmap='jet', alpha=0.6) #scale=120, alpha=0.9, width=0.003
    #Rgdf.plot(ax=ax1, facecolor='lightgrey', edgecolor='lightgrey', alpha=0.2)

    ax1.set_title('A: Wind Velocity Magnitude (m/s)', loc='left', pad=15, weight='bold')

    # --- MAP B: Flow (Quiver) ---
    QUIVsampled_df = ped_df.iloc[::5].copy()
    # Plot binned_df instead of the raw xx rows
    #qv = ax2.quiver(binned_df['x_bin'], binned_df['y_bin'], binned_df['U'], binned_df['V'], binned_df['u_mag'], cmap='jet', scale=120, alpha=0.9, width=0.003) 
    qv = ax2.quiver(QUIVsampled_df['X'], QUIVsampled_df['Y'], QUIVsampled_df['U'], QUIVsampled_df['V'], QUIVsampled_df['u_mag'], cmap='jet', scale=120, alpha=0.9, width=0.003)
    gdf.plot(ax=ax2, facecolor='lightgrey', edgecolor='black', alpha=0.4)
    #Lgdf.plot(ax=ax2, facecolor='green', edgecolor='green', alpha=0.4)

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

    plt.suptitle(f'Pedestrian Wind at 1.5m - {title_suffix}', fontsize=16, y=0.98)
    
    return fig, (ax1, ax2)

def plot_wind_vectors(ped_df, gdf, prk, trees, center_x_utm, center_y_utm, radius=400, title_suffix=""):
    """
    Plots a main wind vector map alongside a separate, dedicated zoom panel 
    for Trafalgar Park with explicit connecting indicator lines.
    Layout: Zoom (Inset) on the LEFT, Main Plot on the RIGHT.
    """
    samp = ped_df.iloc[::5]
    
    # 1. Initialize side-by-side plots (ax1 is Left, ax2 is Right)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={'width_ratios': [1, 1]})

    # Helper function to render spatial layers
    def _plot_layer(target_ax, q_scale=120, q_width=0.003):
        gdf.plot(ax=target_ax, facecolor='lightgrey', edgecolor='black', alpha=0.4)
        target_ax.scatter(trees.geometry.x, trees.geometry.y, marker='o', color='brown', s=15, alpha=0.5, linewidths=2)
        target_ax.scatter(trees.geometry.x, trees.geometry.y, marker='o', color='brown', s=2, alpha=0.5, linewidths=1)


        return target_ax.quiver(
            samp['X'], samp['Y'], samp['U'], samp['V'], samp['u_mag'],
            cmap='jet', scale=q_scale, width=q_width, alpha=0.9
        )

    # 2. Render Main Map on the RIGHT (ax2)
    qv = _plot_layer(ax2)
    ax2.set_xlim(center_x_utm - radius, center_x_utm + radius)
    ax2.set_ylim(center_y_utm - radius, center_y_utm + radius)
    ax2.set_aspect('equal', adjustable='box')
    ax2.set_title(f'Vector Field at  - {title_suffix}', fontsize=12, weight='bold')
    ax2.set_xlabel('Easting (m)')
    ax2.set_ylabel('Northing (m)')
    # Colorbar attached to the main map on the right
    fig.colorbar(qv, ax=ax2, shrink=0.7, label='Wind Speed (m/s)')

    # 3. Render Zoom Panel (Trafalgar Park) on the LEFT (ax1)
    park = prk[prk['name'].str.contains('Trafalgar Park', case=False, na=False)]
    if not park.empty:
        _plot_layer(ax1, q_scale=80, q_width=0.005)
        
        # Calculate bounding box with a tight 30m padding
        minx, miny, maxx, maxy = park.total_bounds
        pad = 30
        xlims, ylims = (minx - pad, maxx + pad), (miny - pad, maxy + pad)
        
        ax1.set_xlim(xlims)
        ax1.set_ylim(ylims)
        ax1.set_aspect('equal', adjustable='box')
        ax1.set_title("Zoom: Trafalgar Park", fontsize=12, weight='bold')
        ax1.set_xlabel('Easting (m)')
        ax1.set_ylabel('Northing (m)')

        # Highlight the target zoom region on the main map (ax2)
        rect = Rectangle((xlims[0], ylims[0]), xlims[1]-xlims[0], ylims[1]-ylims[0], 
                         linewidth=1.5, edgecolor='red', facecolor='none', linestyle='--')
        ax2.add_patch(rect)

        # Draw physical connection lines linking ax1 boundaries to ax2 target region corners
        # Upper connection line (Top-right of ax1 panel to Top-left of ax2 box)
        conn_top = ConnectionPatch(xyA=(xlims[1], ylims[1]), xyB=(xlims[0], ylims[1]), 
                                   coordsA="data", coordsB="data", axesA=ax1, axesB=ax2, 
                                   color="gray", linestyle="--", linewidth=1)
        # Lower connection line (Bottom-right of ax1 panel to Bottom-left of ax2 box)
        conn_bot = ConnectionPatch(xyA=(xlims[1], ylims[0]), xyB=(xlims[0], ylims[0]), 
                                   coordsA="data", coordsB="data", axesA=ax1, axesB=ax2, 
                                   color="gray", linestyle="--", linewidth=1)
        
        fig.add_artist(conn_top)
        fig.add_artist(conn_bot)
    else:
        # Hide the zoom axis cleanly if data isn't found
        ax1.axis('off')
        print("Warning: 'Trafalgar Park' not found in GeoDataFrame.")

    plt.tight_layout()
    return fig, (ax1, ax2)

def get_sun_position(lat: float, lon: float, dt: datetime) -> tuple[float, float]:
    """
    Returns (azimuth_deg, altitude_deg).

    dt must be timezone-aware. UTC offset is derived from dt.utcoffset(),
    so passing Africa/Johannesburg-aware datetimes handles SAST (UTC+2)
    automatically — no hardcoding needed.

    Accuracy: ~0.5° — suitable for urban shadow studies.
    """
    if dt.tzinfo is None:
        raise ValueError("dt must be timezone-aware. "
                         "Use e.g. datetime(..., tzinfo=ZoneInfo('Africa/Johannesburg'))")

    # UTC offset in decimal hours (handles DST automatically)
    utc_offset_h = dt.utcoffset().total_seconds() / 3600.0

    day_of_year = dt.timetuple().tm_yday
    hour_utc    = dt.hour + dt.minute / 60.0 + dt.second / 3600.0
    hour_local  = hour_utc + utc_offset_h          # local clock time

    # Solar declination
    declination = 23.45 * math.sin(math.radians(360 / 365 * (day_of_year - 81)))

    # Equation of Time (minutes)
    b   = math.radians(360 / 364 * (day_of_year - 81))
    eot = 9.87 * math.sin(2 * b) - 7.53 * math.cos(b) - 1.5 * math.sin(b)

    # Local Solar Time
    lstm = 15.0 * utc_offset_h                     # standard meridian for this offset
    tc   = 4.0 * (lon - lstm) + eot                # time correction (minutes)
    lst  = hour_local + tc / 60.0                  # local solar time (hours)

    # Hour Angle
    hra = 15.0 * (lst - 12.0)

    lat_rad = math.radians(lat)
    dec_rad = math.radians(declination)
    hra_rad = math.radians(hra)

    # Altitude
    sin_elev = (math.sin(lat_rad) * math.sin(dec_rad) +
                math.cos(lat_rad) * math.cos(dec_rad) * math.cos(hra_rad))
    sin_elev = max(-1.0, min(1.0, sin_elev))       # clamp for floating-point safety
    elevation = math.asin(sin_elev)
    elevation_deg = math.degrees(elevation)

    # Azimuth — use atan2 to avoid acos domain errors at zenith / horizon
    cos_elev = math.cos(elevation)
    if cos_elev < 1e-10:                           # sun essentially at zenith
        return 0.0, elevation_deg

    cos_az = ((math.sin(dec_rad) * math.cos(lat_rad) -
               math.cos(dec_rad) * math.sin(lat_rad) * math.cos(hra_rad)) / cos_elev)
    cos_az = max(-1.0, min(1.0, cos_az))
    azimuth_deg = math.degrees(math.acos(cos_az))

    if hra > 0:                                    # afternoon — sun west of south
        azimuth_deg = 360.0 - azimuth_deg

    return azimuth_deg, elevation_deg

def calculate_shadows(blds, trees, sun_azimuth, sun_altitude, min_altitude=2.0):
    """
    Returns shadows as a list of (polygon, casting_building_height) tuples
    so downstream code can skip shade for points taller than the caster.
    """
    if sun_altitude <= 0:
        return []

    alt = max(sun_altitude, min_altitude)
    s_factor = 1.0 / math.tan(math.radians(alt))
    angle_rad = math.radians(sun_azimuth + 180)
    dx_unit = math.sin(angle_rad)
    dy_unit = math.cos(angle_rad)

    shadows = []
    for _, bld in blds.iterrows():
        height = bld.get('building_height', 4.0)
        footprint = bld.geometry
        if footprint is None or footprint.is_empty:
            continue
        dx = dx_unit * s_factor * height
        dy = dy_unit * s_factor * height
        shadow_ext = translate(footprint, xoff=dx, yoff=dy)
        full_shadow = unary_union([footprint, shadow_ext])
        shadows.append((full_shadow, height))   # carry the caster height

    for _, bld in trees.iterrows():
        height = bld.get('height')
        footprint = bld.geometry
        if footprint is None or footprint.is_empty:
            continue
        
        #v_type = str(row['type']).lower().strip()
        v_type = str(bld.get('height')).lower().strip()

        # Apply botanical allometric logic using scalar variables
        if 'palm' in v_type:
            c_dia_val = 4.0
        else:
            # Broadleaf standard fallback logic
            c_dia_val = max(3.0, height * 0.55)           # Crown scales with height

        footprint = footprint.buffer(c_dia_val)
        dx = dx_unit * s_factor * height
        dy = dy_unit * s_factor * height
        shadow_ext = translate(footprint, xoff=dx, yoff=dy)
        full_shadow = unary_union([footprint, shadow_ext])
        shadows.append((full_shadow, height))   # carry the caster height

    return shadows

    
def apply_solar_to_df(df: pd.DataFrame, shadows: list) -> pd.DataFrame:
    """
    Height-aware shade classification.
    A point at height Z is only shaded if the casting building is taller than Z.

    shadows: list of (polygon, caster_height) tuples from calculate_shadows_3d.
    """
    df = df.copy()
    df['is_shaded'] = False

    if not shadows:
        return df

    coords = df[['X', 'Y']].to_numpy()
    pts_arr = shapely.points(coords[:, 0], coords[:, 1])
    z_arr = df['Z'].to_numpy()

    # Accumulate shade mask across all buildings
    shade_mask = np.zeros(len(df), dtype=bool)

    for shadow_poly, caster_height in shadows:
        # Only consider points below the caster height
        height_eligible = z_arr < caster_height
        if not height_eligible.any():
            continue
        in_shadow = shapely.contains_properly(shadow_poly, pts_arr)
        shade_mask |= (in_shadow & height_eligible)

    df['is_shaded'] = shade_mask
    return df

def build_utci_layer(summerPed: pd.DataFrame, shadows: list, ta: float, rh: float) -> pd.DataFrame:
    """
    Physics-correct pipeline:
      - Shade is a ground-level question: is the pedestrian standing in shadow?
        → classify shade on the 1.2–1.9m slice (actual pedestrian level)
      - Wind at 10m is a CFD sampling artefact, not a physical height.
        → extract u_mag from the 10m slice and join it down to ground points
      - UTCI is then computed at ground level with 10m wind + ground shade
    """
    # 1. Ground-level shade — this is where the person actually stands
    #ground = ped_df_full[
    #    (ped_df_full['Z'] >= z_ground_min) &
    #    (ped_df_full['Z'] <= z_ground_max)
    #].copy()
    ground = summerPed
    
    #ground = summerPed
    ground = apply_solar_to_df(ground, shadows)   # shade check at Z~1.5m, no height filter needed

    # 2. 10m wind — spatial proxy for pedestrian-level wind exposure
    #wind10m = ped_df_full[
    #    (ped_df_full['Z'] >= z_wind_min) &
    #    (ped_df_full['Z'] <= z_wind_max)
    #][['X', 'Y', 'u_mag']].copy()
    
    wind10m = summerPed#['U10_equivalent']
    #wind10m = summer10m

    # 3. Spatially join wind to ground points (nearest 10m cell → ground point)
    #    Round XY to nearest metre so the merge is tolerant of small offsets
    #    between the two CFD slices
    def _round_xy(df, decimals=0):
        df = df.copy()
        df['_gx'] = np.round(df['X'], decimals)
        df['_gy'] = np.round(df['Y'], decimals)
        return df

    ground  = _round_xy(ground)
    wind10m = _round_xy(wind10m).rename(columns={'U10_equivalent': 'u_mag_10m'})

    merged = ground.merge(
        wind10m[['_gx', '_gy', 'u_mag_10m']],
        on=['_gx', '_gy'],
        how='left',
    ).drop(columns=['_gx', '_gy'])

    # Fall back to ground-level u_mag for any points with no 10m match
    merged['u_mag_10m'] = merged['u_mag_10m'].fillna(merged['u_mag'])

    # 4. UTCI at ground level, using 10m wind
    #def _utci_row(row):
    #    mrt = ta + mrt_offset if not row['is_shaded'] else ta
    #    return calculate_utci_robust(ta, mrt, row['u_mag_10m'], rh)

    #merged['utci'] = merged.apply(_utci_row, axis=1)
    return merged

def plot_utci_summer(summer, buildings_df, trees, center_x_utm, center_y_utm, radius=400, title_suffix=" "):
    """
    Generates a single, high-resolution UTCI Thermal Stress plot for Summer.
    """
    # 1. Determine local min and max for the summer dataset
    vmin = summer['utci'].min()
    vmax = summer['utci'].max()
    
    # 2. Create a clean set of 20 levels
    shared_levels = np.linspace(vmin, vmax, 20)
    
    # 3. Setup a single-plot figure (1 row, 1 col)
    fig, ax1 = plt.subplots(1, 1, figsize=(8, 9))

    # --- Tricontour Contour Fill ---
    cntr = ax1.tricontourf(
        summer['X'], summer['Y'], summer['utci'], 
        levels=shared_levels, cmap='RdYlBu_r', alpha=0.7
    )
    
    # Plot spatial geometries (buildings, terrain layers, etc.)
    #plot_geometries(buildings_df, ax=ax1, facecolor='lightgrey', edgecolor='black', alpha=0.8)
    buildings_df.plot(ax=ax1, facecolor='lightgrey', edgecolor='lightgrey', alpha=0.7)
    #ax1.scatter(trees.geometry.x, trees.geometry.y, marker='*', color='green', s=15, alpha=0.7)#, label='Turbines')
    
    # Optional: If you want to overlay your tree canopy locations onto the axis
    # plot_geometries(trees, ax=ax1, facecolor='none', edgecolor='green', alpha=0.6)

    # --- AXIS & SPATIAL STYLING ---
    ax1.set_aspect('equal', adjustable='box')
    ax1.set_xlabel('Easting (m)')
    ax1.set_ylabel('Northing (m)')
    
    # Force strict AOI clip around your coordinate center
    ax1.set_xlim(center_x_utm - radius, center_x_utm + radius)
    ax1.set_ylim(center_y_utm - radius, center_y_utm + radius)

    # --- COLORBAR STYLING ---
    # Using 'fraction' keeps it nicely scaled to the single plot height
    cbar = fig.colorbar(cntr, ax=ax1, fraction=0.046, pad=0.04)
    cbar.set_label('UTCI (°C) Thermal Stress', weight='bold')

    plt.title(f'Summer UTCI Thermal Stress Mapping {title_suffix}', fontsize=14, pad=15)#, weight='bold')
    plt.tight_layout()
    
    return fig, ax1
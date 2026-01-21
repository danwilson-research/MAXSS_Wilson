# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 14:44:59 2026

@author: dw557
"""

#Code to check for differences between flux engine model runs
#In this case between my FluxEngine output and that by Dan Ford
#this script works with data that has been regridded onto 1deg grid
#code created with assistance from Gemini AI

#Import required packages
import numpy as np
import netCDF4 as nc
import xarray as xr
import matplotlib.pyplot as plt
import os


# --- 1. Load Data ---
dan_f_file = 'E:/MAXSS_working_directory/Ford_et_al_GBC_fco2/flux/2010/06/OceanFluxGHG-month06-jun-2010-v0.nc'

#When using a file that has already been regridded to 1 deg
maxss_file = 'E:/MAXSS_working_directory/output/MAXSS_RUN/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/MAXSS_2010176N16278_AL012010_ALEX_DOY_166_DATE_2010_06_15.nc'

dan_f_ds = nc.Dataset(dan_f_file)
maxss_ds = nc.Dataset(maxss_file)
maxss_xr = xr.open_dataset(maxss_file)

# --- 2. Create the Mask & Spatial Bounds ---
# Create boolean mask from MAXSS 'OF' variable
valid_mask = maxss_xr['OF'][0,:,:].notnull().values

# Get Lat/Lon vectors from MAXSS
maxss_lats = maxss_xr['latitude'].values
maxss_lons = maxss_xr['longitude'].values

# --- CALCULATE ZOOM BOUNDS ---
# Find indices where valid data exists
valid_y_indices, valid_x_indices = np.where(valid_mask)

# Get the actual Lat/Lon values for these edges
# (Using min/max indices to slice the coordinate arrays)
bounds_lat_min = maxss_lats[valid_y_indices.min()]
bounds_lat_max = maxss_lats[valid_y_indices.max()]
bounds_lon_min = maxss_lons[valid_x_indices.min()]
bounds_lon_max = maxss_lons[valid_x_indices.max()]

# Add a buffer (e.g., 5 degrees) to the zoom area
buffer = 5
zoom_lat_min = bounds_lat_min - buffer
zoom_lat_max = bounds_lat_max + buffer
zoom_lon_min = bounds_lon_min - buffer
zoom_lon_max = bounds_lon_max + buffer

# Full grid bounds for slicing Dan F
min_lat, max_lat = maxss_lats.min(), maxss_lats.max()
min_lon, max_lon = maxss_lons.min(), maxss_lons.max()

# --- 3. Calculate Slice Indices for Dan F ---
dan_lats = dan_f_ds['latitude'][:]
dan_lons = dan_f_ds['longitude'][:]

lat_indices = np.where((dan_lats >= min_lat) & (dan_lats <= max_lat))[0]
lon_indices = np.where((dan_lons >= min_lon) & (dan_lons <= max_lon))[0]

lat_start, lat_end = lat_indices[0], lat_indices[-1] + 1
lon_start, lon_end = lon_indices[0], lon_indices[-1] + 1

# Check orientation for flipping logic
is_dan_descending = dan_f_ds['latitude'][0] > dan_f_ds['latitude'][-1]
is_maxss_descending = maxss_lats[0] > maxss_lats[-1]

# --- 4. Loop Through Variables and Plot ---
var_list = maxss_ds.variables
output_dir = 'E:/MAXSS_working_directory/plots/three_panel_zoom/'
os.makedirs(output_dir, exist_ok=True)

for var_name in var_list:
    # Skip if not in Dan F
    if var_name not in dan_f_ds.variables:
        continue

    # Metadata
    long_name = getattr(maxss_ds.variables[var_name], 'long_name', var_name)
    units = getattr(maxss_ds.variables[var_name], 'units', '')
    x_label = f"{long_name} ({units})" if units else long_name

    # --- Extract MAXSS Data ---
    data_maxss = maxss_ds[var_name][:]
    
    # Handle 3D data (Time, Lat, Lon) -> Force to 2D
    if data_maxss.ndim == 3:
        data_maxss = data_maxss[0, :, :]
    
    data_maxss = np.squeeze(data_maxss)
    
    if data_maxss.ndim != 2: 
        continue
    
    data_maxss_flat = data_maxss.flatten()
    data_maxss_clean = data_maxss_flat[~np.isnan(data_maxss_flat)]

    # --- Extract Dan F Data ---
    try:
        if dan_f_ds[var_name].ndim == 3:
            data_dan = dan_f_ds[var_name][:, lat_start:lat_end, lon_start:lon_end]
            if data_dan.shape[0] == 1:
                data_dan = data_dan[0, :, :]
            else:
                data_dan = data_dan[0, :, :]
        elif dan_f_ds[var_name].ndim == 2:
            data_dan = dan_f_ds[var_name][lat_start:lat_end, lon_start:lon_end]
        else:
            continue
            
        data_dan = np.squeeze(data_dan)

        if is_dan_descending != is_maxss_descending:
             data_dan = np.flipud(data_dan)

        data_dan_masked = np.where(valid_mask, data_dan, np.nan)
        data_dan_flat = data_dan_masked.flatten()
        data_dan_clean = data_dan_flat[~np.isnan(data_dan_flat)]
        
    except Exception as e:
        print(f"Error processing {var_name}: {e}")
        continue

    # --- Print Counts ---
    count_maxss = len(data_maxss_clean)
    count_dan = len(data_dan_clean)
    print(f"Variable: {var_name: <20} | Valid Counts -> MAXSS: {count_maxss}, Dan F: {count_dan}")

    if count_maxss == 0 and count_dan == 0:
        continue

    # --- Plotting (3 Panels) ---
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Determine color scale
    combined_clean = np.concatenate([data_maxss_clean, data_dan_clean])
    if len(combined_clean) > 0:
        vmin, vmax = np.min(combined_clean), np.max(combined_clean)
    else:
        vmin, vmax = None, None

    # Panel 1: MAXSS Map
    im1 = axes[0].pcolormesh(maxss_lons, maxss_lats, data_maxss, 
                             cmap='viridis', vmin=vmin, vmax=vmax, shading='auto')
    axes[0].set_title(f"MAXSS: {var_name}")
    axes[0].set_xlabel("Longitude")
    axes[0].set_ylabel("Latitude")
    axes[0].set_xlim(zoom_lon_min, zoom_lon_max) # APPLY ZOOM
    axes[0].set_ylim(zoom_lat_min, zoom_lat_max) # APPLY ZOOM
    axes[0].grid(True, linestyle='--', alpha=0.5, color='gray') # ADD GRID
    plt.colorbar(im1, ax=axes[0], orientation='vertical', fraction=0.046, pad=0.04)

    # Panel 2: Dan F Map
    im2 = axes[1].pcolormesh(maxss_lons, maxss_lats, data_dan_masked, 
                             cmap='viridis', vmin=vmin, vmax=vmax, shading='auto')
    axes[1].set_title(f"Dan F (Masked): {var_name}")
    axes[1].set_xlabel("Longitude")
    axes[1].set_xlim(zoom_lon_min, zoom_lon_max) # APPLY ZOOM
    axes[1].set_ylim(zoom_lat_min, zoom_lat_max) # APPLY ZOOM
    axes[1].grid(True, linestyle='--', alpha=0.5, color='gray') # ADD GRID
    plt.colorbar(im2, ax=axes[1], orientation='vertical', fraction=0.046, pad=0.04)

    # Panel 3: Histogram (Density)
    if vmin is not None and vmax is not None and vmin != vmax:
        bins = np.linspace(vmin, vmax, 50)
    else:
        bins = 50
        
    # Using density=True
    axes[2].hist(data_maxss_clean, bins=bins, alpha=0.5, label=f'MAXSS', 
                 color='blue', density=True, edgecolor='black', linewidth=0.5)
    axes[2].hist(data_dan_clean, bins=bins, alpha=0.5, label=f'Dan F', 
                 color='orange', density=True, edgecolor='black', linewidth=0.5)
    
    axes[2].set_title(f"Histogram: {long_name}")
    axes[2].set_xlabel(x_label)
    axes[2].set_ylabel("Density")
    axes[2].legend()
    axes[2].grid(True, alpha=0.3, linestyle='--')

    plt.tight_layout()
    
    safe_var_name = var_name.replace(" ", "_").replace("/", "-")
    filename = f"{output_dir}panel_comparison_{safe_var_name}.png"
    plt.savefig(filename)
    plt.close()

print("Processing complete. Plots saved.")
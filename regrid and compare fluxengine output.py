# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 17:40:32 2026

@author: dw557
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os

# --- 1. Load Data ---
dan_f_file = 'E:/MAXSS_working_directory/Ford_et_al_GBC_fco2/flux/2010/06/OceanFluxGHG-month06-jun-2010-v0.nc'
maxss_file = 'E:/MAXSS_working_directory/output/MAXSS_RUN/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/MAXSS_2010176N16278_AL012010_ALEX_DOY_166_DATE_2010_06_15.nc'

# Load both as Xarray Datasets immediately
ds_dan = xr.open_dataset(dan_f_file)
ds_maxss = xr.open_dataset(maxss_file)

# --- 2. Normalize Longitudes (Crucial for Regridding) ---
# This ensures both datasets use -180 to 180 range before we try to match them.
def normalize_longitude(ds):
    # Adjusts 0..360 to -180..180
    ds.coords['longitude'] = (ds.coords['longitude'] + 180) % 360 - 180
    ds = ds.sortby(ds.longitude)
    return ds

ds_dan = normalize_longitude(ds_dan)
ds_maxss = normalize_longitude(ds_maxss)

# --- 3. REGRID MAXSS TO DAN F GRID ---
print("Regridding MAXSS data onto Dan F grid...")
# This is the magic line. It interpolates MAXSS vars to the Lat/Lon points of Dan F.
# method='nearest' preserves integer flags; 'linear' is better for smooth fields like temp/pressure.
# We use 'nearest' here to be safe with masks/boundaries, but you can swap to 'linear'.
ds_maxss_regridded = ds_maxss.interp(
    latitude=ds_dan.latitude, 
    longitude=ds_dan.longitude, 
    method='nearest' 
)

# --- 4. Define Zoom Area (Based on where MAXSS data actually exists) ---
# We look for valid data in the 'OF' variable (or the first variable available)
valid_mask = ds_maxss_regridded['OF'].notnull()
if 'time' in valid_mask.dims:
    valid_mask = valid_mask.isel(time=0)

# Get lat/lon arrays of the target grid (Dan F)
grid_lats = ds_dan['latitude'].values
grid_lons = ds_dan['longitude'].values

# Find the bounding box of the valid data on this new grid
valid_y_indices, valid_x_indices = np.where(valid_mask.values)

if len(valid_y_indices) > 0:
    # Calculate bounds with a buffer
    buffer = 5
    zoom_lat_min = grid_lats[valid_y_indices.min()] - buffer
    zoom_lat_max = grid_lats[valid_y_indices.max()] + buffer
    zoom_lon_min = grid_lons[valid_x_indices.min()] - buffer
    zoom_lon_max = grid_lons[valid_x_indices.max()] + buffer
else:
    # Fallback if mask is empty
    zoom_lat_min, zoom_lat_max = -90, 90
    zoom_lon_min, zoom_lon_max = -180, 180

# --- 5. Loop and Plot ---
var_list = ds_maxss.variables
output_dir = 'E:/MAXSS_working_directory/plots/regridded_comparison/'
os.makedirs(output_dir, exist_ok=True)

for var_name in var_list:
    # Skip if not in Dan F or if it's a coordinate
    if var_name not in ds_dan.variables or var_name in ['latitude', 'longitude', 'time']:
        continue

    # Metadata
    long_name = getattr(ds_maxss[var_name], 'long_name', var_name)
    units = getattr(ds_maxss[var_name], 'units', '')
    x_label = f"{long_name} ({units})" if units else long_name

    # --- Extract Data (Both are now on Dan F grid) ---
    # 1. MAXSS (Regridded)
    data_maxss = ds_maxss_regridded[var_name]
    if 'time' in data_maxss.dims: data_maxss = data_maxss.isel(time=0)
    data_maxss = data_maxss.values # Convert to numpy

    # 2. Dan F (Original Grid)
    data_dan = ds_dan[var_name]
    if 'time' in data_dan.dims: data_dan = data_dan.isel(time=0)
    data_dan = data_dan.values # Convert to numpy

    # --- Apply Masking ---
    # We want to focus on the area where MAXSS has data.
    # We create a combined mask: Valid in MAXSS AND Valid in Dan F?
    # Or just mask by MAXSS extent (so we can see if Dan F covers it).
    
    # Mask Dan F data outside the MAXSS study area (using the mask we built earlier)
    # Note: data_dan is full global, we mask it to the 'valid_mask' region for visual comparison
    data_dan_masked = np.where(valid_mask.values, data_dan, np.nan)
    data_maxss_masked = np.where(valid_mask.values, data_maxss, np.nan)

    # Flatten for Histogram
    flat_maxss = data_maxss_masked.flatten()
    clean_maxss = flat_maxss[~np.isnan(flat_maxss)]

    flat_dan = data_dan_masked.flatten()
    clean_dan = flat_dan[~np.isnan(flat_dan)]

    # --- Print Counts ---
    print(f"Variable: {var_name: <20} | Valid Pixels (Intersection) -> MAXSS: {len(clean_maxss)}, Dan F: {len(clean_dan)}")

    if len(clean_maxss) == 0 and len(clean_dan) == 0:
        continue

    # --- Plotting ---
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Determine Color Scale
    combined_clean = np.concatenate([clean_maxss, clean_dan])
    if len(combined_clean) > 0:
        vmin, vmax = np.min(combined_clean), np.max(combined_clean)
    else:
        vmin, vmax = None, None

    # Panel 1: MAXSS (Regridded)
    im1 = axes[0].pcolormesh(grid_lons, grid_lats, data_maxss_masked, 
                             cmap='viridis', vmin=vmin, vmax=vmax, shading='auto')
    axes[0].set_title(f"MAXSS (Regridded to Dan F): {var_name}")
    axes[0].set_xlabel("Longitude")
    axes[0].set_ylabel("Latitude")
    axes[0].set_xlim(zoom_lon_min, zoom_lon_max)
    axes[0].set_ylim(zoom_lat_min, zoom_lat_max)
    axes[0].grid(True, linestyle='--', alpha=0.5)
    plt.colorbar(im1, ax=axes[0])

    # Panel 2: Dan F (Original)
    im2 = axes[1].pcolormesh(grid_lons, grid_lats, data_dan_masked, 
                             cmap='viridis', vmin=vmin, vmax=vmax, shading='auto')
    axes[1].set_title(f"Dan F (Masked to Area): {var_name}")
    axes[1].set_xlabel("Longitude")
    axes[1].set_xlim(zoom_lon_min, zoom_lon_max)
    axes[1].set_ylim(zoom_lat_min, zoom_lat_max)
    axes[1].grid(True, linestyle='--', alpha=0.5)
    plt.colorbar(im2, ax=axes[1])

    # Panel 3: Histogram
    bins = np.linspace(vmin, vmax, 50) if (vmin is not None and vmax is not None) else 50
    
    axes[2].hist(clean_maxss, bins=bins, alpha=0.5, label='MAXSS', 
                 color='blue', density=True, edgecolor='black', linewidth=0.5)
    axes[2].hist(clean_dan, bins=bins, alpha=0.5, label='Dan F', 
                 color='orange', density=True, edgecolor='black', linewidth=0.5)
    
    axes[2].set_title(f"Density Histogram: {long_name}")
    axes[2].set_xlabel(x_label)
    axes[2].legend()
    axes[2].grid(True, alpha=0.3, linestyle='--')

    plt.tight_layout()
    safe_var_name = var_name.replace(" ", "_").replace("/", "-")
    plt.savefig(f"{output_dir}regrid_compare_{safe_var_name}.png")
    plt.close()

print("Processing complete. Regridded plots saved.")
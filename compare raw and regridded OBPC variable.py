# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 14:38:14 2026

@author: dw557
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# 1. Load the Datasets
raw_input_df = "E:/MAXSS_working_directory/Ford_et_al_GBC_fco2/flux/2010/06/OceanFluxGHG-month06-jun-2010-v0.nc"
resampled_output_df = "E:/MAXSS_working_directory/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/Resampled_for_fluxengine_Ford_et_al_pco2.nc"

try:
    raw_input = nc.Dataset(raw_input_df)
    resampled_output = nc.Dataset(resampled_output_df)

    # 2. Load the Variables
    # Using [0, :, :] assumes the dimensions are (time, lat, lon) or similar
    raw_data = raw_input['OBPC'][0, :, :]
    resampled_data = resampled_output['pCO2water_mean'][0, :, :]

    # 3. Load Coordinates (Lat/Lon)
    # Note: Variable names for coordinates might vary (e.g., 'lat', 'latitude'). 
    # Adjust 'lat'/'lon' below if your file uses different names.
    
    # Check for raw input coordinates
    if 'lat' in raw_input.variables:
        raw_lat = raw_input['lat'][:]
        raw_lon = raw_input['lon'][:]
    else:
        raw_lat = raw_input['latitude'][:]
        raw_lon = raw_input['longitude'][:]

    # Check for resampled output coordinates
    if 'lat' in resampled_output.variables:
        res_lat = resampled_output['lat'][:]
        res_lon = resampled_output['lon'][:]
    else:
        res_lat = resampled_output['latitude'][:]
        res_lon = resampled_output['longitude'][:]

    # 4. Determine Geographical Extent from Resampled Output
    # Extent format: [min_lon, max_lon, min_lat, max_lat]
    extent = [np.min(res_lon), np.max(res_lon), np.min(res_lat), np.max(res_lat)]

    # 5. Determine Common Color Scale (vmin, vmax)
    # We use nanmin/nanmax to ignore missing values
    vmin = min(np.nanmin(raw_data), np.nanmin(resampled_data))
    vmax = max(np.nanmax(raw_data), np.nanmax(resampled_data))

    # 6. Plotting
    fig = plt.figure(figsize=(14, 6))

    # -- Plot 1: Raw Input --
    ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
    ax1.set_extent(extent, crs=ccrs.PlateCarree())
    
    # Add map features
    ax1.add_feature(cfeature.COASTLINE)
    ax1.add_feature(cfeature.BORDERS, linestyle=':')
    ax1.gridlines(draw_labels=True)

    # Plot data
    # pcolormesh handles 1D lat/lon arrays automatically matching the 2D data dimensions
    im1 = ax1.pcolormesh(raw_lon, raw_lat, raw_data, 
                         transform=ccrs.PlateCarree(), 
                         vmin=vmin, vmax=vmax, cmap='viridis', shading='auto')
    ax1.set_title('Raw Input (OBPC)')

    # -- Plot 2: Resampled Output --
    ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())
    ax2.set_extent(extent, crs=ccrs.PlateCarree())
    
    # Add map features
    ax2.add_feature(cfeature.COASTLINE)
    ax2.add_feature(cfeature.BORDERS, linestyle=':')
    ax2.gridlines(draw_labels=True)

    # Plot data
    im2 = ax2.pcolormesh(res_lon, res_lat, resampled_data, 
                         transform=ccrs.PlateCarree(), 
                         vmin=vmin, vmax=vmax, cmap='viridis', shading='auto')
    ax2.set_title('Resampled Output (pCO2water_mean)')

    # -- Shared Colorbar --
    # Add a colorbar on the right side of the figure
    cbar = fig.colorbar(im2, ax=[ax1, ax2], orientation='vertical', fraction=0.046, pad=0.04)
    cbar.set_label('pCO2')

    plt.suptitle("Comparison of Raw vs Resampled Data")
    plt.show()

    # Close datasets
    raw_input.close()
    resampled_output.close()

except FileNotFoundError:
    print("Error: Files not found. Please check the file paths.")
except Exception as e:
    print(f"An error occurred: {e}")


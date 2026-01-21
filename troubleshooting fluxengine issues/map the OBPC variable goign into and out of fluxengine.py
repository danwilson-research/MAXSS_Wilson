# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 10:47:45 2026

@author: dw557
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# 1. Load the Datasets
fluxengine_input_df = "E:/MAXSS_working_directory/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/Resampled_for_fluxengine_Ford_et_al_pco2.nc"
fluxengine_output_df = "E:/MAXSS_working_directory/output/MAXSS_RUN/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/MAXSS_2010176N16278_AL012010_ALEX_DOY_166_DATE_2010_06_15.nc"

# Load in the netcdf files
fluxengine_input = nc.Dataset(fluxengine_input_df)
fluxengine_output = nc.Dataset(fluxengine_output_df)

# Extract Data variables
# Using np.squeeze to remove single-dimensional entries (like the time dimension) if necessary
fe_input_obpc = np.squeeze(fluxengine_input['pCO2water_mean'][0,:,:])
fe_output_obpc = np.squeeze(fluxengine_output['OBPC'][0,:,:])

# 2. Extract Coordinates (Latitude and Longitude)
# Assuming standard variable names. Adjust 'lat'/'lon' if your file uses 'latitude'/'longitude'
try:
    lats = fluxengine_input.variables['lat'][:]
    lons = fluxengine_input.variables['lon'][:]
except KeyError:
    try:
        lats = fluxengine_input.variables['latitude'][:]
        lons = fluxengine_input.variables['longitude'][:]
    except KeyError:
        print("Error: Could not find latitude/longitude variables. Please check variable names.")
        raise

# Handle 1D vs 2D coordinates for meshgrid
if lats.ndim == 1 and lons.ndim == 1:
    lons_grid, lats_grid = np.meshgrid(lons, lats)
else:
    lons_grid, lats_grid = lons, lats

# 3. Setup Plotting
# Determine common color scale (vmin, vmax) for direct comparison
vmin = min(np.nanmin(fe_input_obpc), np.nanmin(fe_output_obpc))
vmax = max(np.nanmax(fe_input_obpc), np.nanmax(fe_output_obpc))

# Create figure with Cartopy Projection
fig = plt.figure(figsize=(14, 6))
projection = ccrs.PlateCarree()

# --- Subplot 1: Input Data ---
ax1 = fig.add_subplot(1, 2, 1, projection=projection)
ax1.coastlines()
ax1.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

# Plot data
mesh1 = ax1.pcolormesh(lons_grid, lats_grid, fe_input_obpc, 
                       transform=projection, cmap='viridis', vmin=vmin, vmax=vmax)

ax1.set_title("Input: pCO2water_mean")

# Gridlines
gl1 = ax1.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl1.top_labels = False
gl1.right_labels = False

# --- Subplot 2: Output Data ---
ax2 = fig.add_subplot(1, 2, 2, projection=projection)
ax2.coastlines()
ax2.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

# Plot data
mesh2 = ax2.pcolormesh(lons_grid, lats_grid, fe_output_obpc, 
                       transform=projection, cmap='viridis', vmin=vmin, vmax=vmax)

ax2.set_title("Output: OBPC")

# Gridlines
gl2 = ax2.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl2.top_labels = False
gl2.right_labels = False
gl2.left_labels = False # Remove left labels for the second plot to reduce clutter

# 4. Add a shared Colorbar
# Adjust position as needed
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7]) # [left, bottom, width, height]
fig.colorbar(mesh2, cax=cbar_ax, label='pCO2 / OBPC')

plt.subplots_adjust(wspace=0.05, right=0.9)
plt.show()

# Close datasets
fluxengine_input.close()
fluxengine_output.close()
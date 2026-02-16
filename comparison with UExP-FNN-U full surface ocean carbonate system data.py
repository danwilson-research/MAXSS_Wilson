# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 12:09:59 2026

@author: dw557
"""

# Script to compare air-sea co2 flux between Dan Wilson FluxENgine output and 
# 'UExP-FNN-U full surface ocean carbonate system' data

# DO NOT FORGET TO ADD IN ICE COVER INTO EQUATION

#Import required packages
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os

#Path to the UExP-FNN-U global data (Use the one large output from dan rather than flux folders)
UEXP_FNN_path = "E:/MAXSS_working_directory/Ford_et_al_GBC_fco2/Fordetal_UExP-FNN-U_surface-carbonate-system_v2025-1.nc"

#Path to directory containing FluxEngine output
fluxengine_dir = "E:/MAXSS_working_directory/output/MAXSS_RUN/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX"

## Extract 2D common land mask ##

#Load in UExP-FNN-U global data
UExP_FNN_raw = xr.open_dataset(UEXP_FNN_path)

## load in example fluxengine data
fluxengine_example_raw = xr.open_dataset('E:/MAXSS_working_directory/output/MAXSS_RUN/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/MAXSS_2010176N16278_AL012010_ALEX_DOY_160_DATE_2010_06_09.nc')

#Transpose UExP to match the (time, lat, lon) order of FluxEngine output
UExP_FNN_example_flux = UExP_FNN_raw.flux.isel(time=0).transpose("latitude", "longitude")
fluxengine_example_flux = fluxengine_example_raw.OF.isel(time=0)

#Perform regridding via nearest neighbour to upscale 1deg UExP_FNN data to 0.25deg
UExP_FNN_example_flux_regridded = UExP_FNN_example_flux.interp_like(fluxengine_example_flux, 
                                                                    method="nearest")

#Now find where valid data (no NaN's located)
UExP_FNN_mask = UExP_FNN_example_flux_regridded.notnull()
fluxengine_mask = fluxengine_example_flux.notnull()

#Create the common mask between the two datasets
common_mask = UExP_FNN_mask & fluxengine_mask.values

###############################################################################

## Sanity check the combined land mask usign a plot

# --- Step 1: Create the Classification Map  ---
# 0 = No Data (Background)
# 1 = FluxEngine Only (Lost UExP)
# 2 = Common Data (Kept)
# -1 = UExP Only (Lost FluxEngine - Rare but possible)

mask_uexp_int = UExP_FNN_mask.astype(int)
mask_fe_int = fluxengine_mask.astype(int)

classification_map = xr.full_like(mask_uexp_int, 0) # Start with 0
# Logic: 
# If FluxEngine has data (1) AND UExP is missing (0) -> Label 1 (Lost UExP)
classification_map = xr.where((mask_fe_int == 1) & (mask_uexp_int == 0), 1, classification_map)
# If FluxEngine is missing (0) AND UExP has data (1) -> Label -1 (Lost FluxEngine)
classification_map = xr.where((mask_fe_int == 0) & (mask_uexp_int == 1), -1, classification_map)
# If BOTH have data -> Label 2 (Kept)
classification_map = xr.where((mask_fe_int == 1) & (mask_uexp_int == 1), 2, classification_map)

# --- Step 2: Plotting with Cartopy ---

# Create figure with PlateCarree projection
fig = plt.figure(figsize=(12, 10))
ax = plt.axes(projection=ccrs.PlateCarree())

# Add Map Features
ax.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='lightgray')
ax.coastlines(resolution='50m', zorder=2) # Higher res coastlines for 0.25 deg data

# Add Gridlines
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

# Define Colormap
# -1: Orange (Lost FluxEngine)
#  0: White/Transparent (No Data)
#  1: Red (Lost UExP - The critical loss)
#  2: Blue (Kept Data)
cmap = mcolors.ListedColormap(['orange', 'none', 'red', 'blue'])
bounds = [-1.5, -0.5, 0.5, 1.5, 2.5]
norm = mcolors.BoundaryNorm(bounds, cmap.N)

# Plot Data
# Note: transform=ccrs.PlateCarree() is crucial so cartopy knows the data is lat/lon
im = classification_map.plot(ax=ax, 
                             cmap=cmap, 
                             norm=norm, 
                             add_colorbar=False, 
                             transform=ccrs.PlateCarree())

#  Custom Legend
# We use a distinct colorbar to act as a legend
cbar = plt.colorbar(im, ticks=[-1, 0, 1, 2], orientation='horizontal', pad=0.05, shrink=0.8)
cbar.ax.set_xticklabels([
    'UExP Only\n(Lost FE)', 
    'No Data', 
    'FluxEngine Only\n(Lost due to coarse UExP)', 
    'Common Data\n(Kept)'
])

ax.set_title("Data Availability: 1° UExP vs 0.25° FluxEngine", fontsize=14)
plt.show()

###############################################################################

## Calculate the area lost through comparing the two datasets

# --- Step 1: Calculate the Area of Each Grid Cell ---
# Earth's radius in km
R = 6371.0 

# Get latitude and longitude resolution (0.25 assumed)
# Infer from the data coordinates to be safe
dlat = np.abs(fluxengine_example_flux.latitude[1] - fluxengine_example_flux.latitude[0])
dlon = np.abs(fluxengine_example_flux.longitude[1] - fluxengine_example_flux.longitude[0])

# Convert degrees to radians
dlon_rad = np.deg2rad(dlon)
dlat_rad = np.deg2rad(dlat)

# Get the latitude of every grid point
# We rely on the FluxEngine grid (0.25 deg) since that is our base resolution
latitudes = fluxengine_example_flux.latitude

# Calculate area for each latitude band
# Formula: Area = R^2 * dlon * (sin(lat_top) - sin(lat_bottom))
# Calculate the edges of the cells (lat +/- half resolution)
lat_top_rad = np.deg2rad(latitudes + dlat/2)
lat_bot_rad = np.deg2rad(latitudes - dlat/2)

# This creates a 1D array of areas (one value per latitude)
cell_areas_1d = (R**2) * dlon_rad * (np.sin(lat_top_rad) - np.sin(lat_bot_rad))

# Broadcast 1D latitude array to the full 2D grid shape (lat, lon)
# This gives us a map where every pixel value is its area in km^2
pixel_area_map, _ = xr.broadcast(cell_areas_1d, fluxengine_example_flux)

# --- Step 2: Sum Areas Based on Your Classification ---

# Total Potential FluxEngine Area (Mask = 1 or 2 in your map)
#    This includes the Red (Lost) and Blue (Kept) pixels
total_fe_mask = (classification_map == 1) | (classification_map == 2)
total_fe_area = pixel_area_map.where(total_fe_mask).sum()

# Area Lost (Mask = 1)
#    This is the Red pixels (FluxEngine valid, but UExP missing)
lost_area_mask = (classification_map == 1)
lost_area = pixel_area_map.where(lost_area_mask).sum()

# Area Kept (Mask = 2)
#    This is the Blue pixels (Common data)
kept_area_mask = (classification_map == 2)
kept_area = pixel_area_map.where(kept_area_mask).sum()

# --- Step 3: Output the Results ---
print("--- Area Statistics ---")
print(f"Total Valid FluxEngine Area:  {total_fe_area.values/1e6:.2f} million km^2")
print(f"Area Retained (Common Data):  {kept_area.values/1e6:.2f} million km^2")
print(f"Area Lost (Coarse Mask):      {lost_area.values/1e6:.2f} million km^2")
print("-" * 30)
print(f"Percentage of Area LOST:      {(lost_area / total_fe_area * 100).values:.2f}%")

###############################################################################


#Get a list of all the fluxengine output files to loop through
fluxengine_files = os.listdir(fluxengine_dir)

print(f"Starting hourly integration for {len(fluxengine_files)} files...")

# Create an empty list to store daily results
hourly_fluxengine_flux = []

#Now begin the process of looping through each timestep
for file in fluxengine_files:
    # Construct full path and open dataset
    file_path = os.path.join(fluxengine_dir, file)
    
    # Using 'with' ensures the file is closed after reading
    with xr.open_dataset(file_path) as fluxengine_hourly:
        # 1. Extract the Flux variable
        fe_flux = fluxengine_hourly.OF
        
        # 2. Apply the Common Mask
        # Sets all pixels outside common area to NaN
        fe_masked = fe_flux.where(common_mask)
        
        # 3. Spatially sum Flux (Result is g C m-2 day-1 * km^2) 
        spatial_integrated_flux = (fe_masked * pixel_area_map).sum(dim=['latitude', 'longitude'])
        
        # 2. Unit conversion to Tg C hr-1
        # g -> Tg is 1e-12
        # km2 -> m2 is 1e6
        # day-1 -> hr-1 is / 24
        hourly_Tg = (spatial_integrated_flux * 1e6 * 1e-12) / 24.0
        
        #Append to list
        hourly_fluxengine_flux.append(hourly_Tg)

# Combine into one timeline
total_fluxengine_hourly_Tg = xr.concat(hourly_fluxengine_flux, dim='time')

plt.plot(total_fluxengine_hourly_Tg)

### now move on to DOING SAME CALCULATAION FOR THE uep_fnn data

#CUT OFF END OF DATA AFTER END OF MODEL RUN
#DONT FORGET INTEGRATNIG SEA ICE/ PROPORTION ICE IMPACT
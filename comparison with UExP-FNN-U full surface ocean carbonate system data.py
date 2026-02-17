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
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
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
##Plot the hourly flux data from the 'common area'

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
        
        # 3. Spatially sum Flux (Result is g C m-2 day-1 * km^2) if less than 1 valid non-nan then value = nan 
            # This helps ID end of FluxEngine model run
        spatial_integrated_flux = (fe_masked * pixel_area_map).sum(dim=['latitude', 'longitude'], min_count=1)
        
        # 2. Unit conversion to Tg C hr-1
        # g -> Tg is 1e-12
        # km2 -> m2 is 1e6
        # day-1 -> hr-1 is / 24
        hourly_Tg = (spatial_integrated_flux * 1e6 * 1e-12) / 24.0
        
        #Append to list
        hourly_fluxengine_flux.append(hourly_Tg)

# Combine into one timeline
fluxengine_hourly_Tg = xr.concat(hourly_fluxengine_flux, dim='time')

#Trim data to remove timesteps where all data = NaN as after end of model run or at the start

#Identify valid data (where Flux is NOT NaN)
valid_mask = fluxengine_hourly_Tg.notnull()

if valid_mask.any():
    # Get all times that have valid data
    valid_times = fluxengine_hourly_Tg.time[valid_mask]
    
    # 2. Find the FIRST and LAST valid timestamps
    first_valid_time = valid_times.min().values
    last_valid_time = valid_times.max().values
    
    print('trimming FluxEngine data to model runtime')
    print(f"  Start: {str(first_valid_time)[:16]}")
    print(f"  End:   {str(last_valid_time)[:16]}")

    # 3. Slice BOTH datasets to this exact window
    # This removes leading NaNs AND trailing NaNs
    fluxengine_hourly_Tg_trimmed = fluxengine_hourly_Tg.sel(time=slice(first_valid_time, last_valid_time))
    

else:
    print("ERROR: No valid data found in the entire time series.")
    fluxengine_hourly_Tg_trimmed = fluxengine_hourly_Tg
    
plt.plot(fluxengine_hourly_Tg_trimmed)
plt.title('FluxEngine hourly flux')
plt.xlabel('timestep')
plt.ylabel('Hourly Flux (Tg C hr-1)')
plt.show()

###############################################################################
### Now move on to calculating hourly flux for UExP 

# 1. Regrid the FULL UExP flux time series to the 0.25deg regional grid
# Using nearest neighbor to match the "blocky" style of the mask
UExP_regridded_full = UExP_FNN_raw.flux.transpose("time", "latitude", "longitude").interp_like(
    fluxengine_example_flux, method="nearest")

# 2. Spatially integrate the FULL UExP timeline first (keep the units consistent)
UExP_masked_full = UExP_regridded_full.where(common_mask)
UExP_spatial_sum = (UExP_masked_full * pixel_area_map).sum(dim=['latitude', 'longitude'])

# Convert UExP (g C m-2 day-1) to (Tg C hr-1)
UExP_Tg_hourly_rate = (UExP_spatial_sum * 1e6 * 1e-12) / 24.0

# 3. Create an EXPLICIT lookup dictionary: {(year, month): value}
# This ensures we only ever pull June 2010 for a June 2010 hour.
uexp_lookup = {
    (int(t.dt.year), int(t.dt.month)): val.item() 
    for t, val in zip(UExP_Tg_hourly_rate.time, UExP_Tg_hourly_rate)}

# 4. Broadcast the baseline to match your storm's hourly timeline
# We use a list comprehension to look up the (year, month) for every hour in TC Alex
uexp_baseline_values = []
for t in fluxengine_hourly_Tg_trimmed.time:
    key = (int(t.dt.year), int(t.dt.month))
    # This will pull the exact matching month or crash if the month is missing (safer!)
    uexp_baseline_values.append(uexp_lookup[key])

# 5. Turn it back into an Xarray DataArray 
uexp_baseline_hourly = xr.DataArray(
    uexp_baseline_values, 
    coords={'time': fluxengine_hourly_Tg_trimmed.time}, 
    dims=['time'])

###############################################################################
## Plot the analysis

# Calculate Total Mass (Sum of hourly rates) ---
# Since units are Tg/hr, the sum of all hours gives the total Tg for the period
fe_total_mass = fluxengine_hourly_Tg_trimmed.sum().values
uexp_total_mass = uexp_baseline_hourly.sum().values
net_difference = fe_total_mass - uexp_total_mass

# Create the Plot 
fig, ax = plt.subplots(figsize=(12, 6))

# Add the Zero Line (Black Dashed)
# zorder=1 ensures it stays behind the data lines
ax.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.8, zorder=1)

# Plot UExP Baseline
ax.plot(uexp_baseline_hourly.time, uexp_baseline_hourly, 
        label='UExP Baseline (Monthly Mean)', color='red', linestyle='--', linewidth=1.5)

# Plot FluxEngine
ax.plot(fluxengine_hourly_Tg_trimmed.time, fluxengine_hourly_Tg_trimmed, 
        label='FluxEngine (TC Alex Hourly)', color='blue', linewidth=1.5)

# Shade the anomaly area
ax.fill_between(fluxengine_hourly_Tg_trimmed.time, uexp_baseline_hourly, fluxengine_hourly_Tg_trimmed, 
                 where=(fluxengine_hourly_Tg_trimmed >= uexp_baseline_hourly),
                 facecolor='green', alpha=0.15, interpolate=True)
ax.fill_between(fluxengine_hourly_Tg_trimmed.time, uexp_baseline_hourly, fluxengine_hourly_Tg_trimmed, 
                 where=(fluxengine_hourly_Tg_trimmed < uexp_baseline_hourly),
                 facecolor='orange', alpha=0.15, interpolate=True)

# --- Step 3: Add the Summary Box in Top Right ---
stats_text = (
    f"Total UExP C flux: {uexp_total_mass:.4f} Tg C\n"
    f"Total MAXSS FluxEngine C flux: {fe_total_mass:.4f} Tg C\n"
    f"---------------------------\n"
    f"MAXSS FluxEngine - UExP = {net_difference:.4f} Tg C"
)

# transform=ax.transAxes ensures (0.95, 0.95) is relative to the plot box
ax.text(0.97, 0.95, stats_text, transform=ax.transAxes, 
        fontsize=10, verticalalignment='top', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))

# --- Step 4: Formatting ---
ax.set_title("Air-Sea CO2 Flux Comparison: TC Alex", fontsize=14, fontweight='bold')
ax.set_ylabel("Flux Rate (Tg C hr$^{-1}$)", fontsize=12)
ax.set_xlabel("Date and Time (UTC)", fontsize=12)

# Format the dates
ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b\n%H:%M'))
plt.xticks(rotation=0)

# 2. Force Y-Axis to Include 0
# Get current limits
ymin, ymax = ax.get_ylim()
# Expand limits if 0 is not already included
if ymin > 0: ax.set_ylim(bottom=0)
if ymax < 0: ax.set_ylim(top=0)

ax.legend(loc='upper left')
ax.grid(True, alpha=0.2)

plt.tight_layout()
plt.show()


#DONT FORGET INTEGRATNIG SEA ICE/ PROPORTION ICE IMPACT
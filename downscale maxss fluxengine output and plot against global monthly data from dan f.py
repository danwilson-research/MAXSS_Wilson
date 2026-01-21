# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 17:06:09 2026

@author: dw557
"""

# Script to make the FluxEngine MAXSS output the same 1 degree resolution as the
# data from Dan Ford

#Import required packages
import os
import netCDF4 as nc
import xarray as xr
import numpy as np

#Load in global monthly flux data
global_monthly_subset = nc.Dataset('E:/MAXSS_working_directory/flux_subset_extracted.nc')

lon_grid = global_monthly_subset['lon'][:]
lat_grid = global_monthly_subset['lat'][:]

# Source directory
fluxengine_data_path = "E:/MAXSS_working_directory/output/MAXSS_RUN/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/"

# New Output directory
# We use r"..." to handle the backslashes safely
output_path = r"E:\MAXSS_working_directory\output\MAXSS_RUN\maxss\storm-atlas\ibtracs\north-atlantic\2010\2010176N16278_AL012010_ALEX\1_degree_resolution_output"

# Create the directory if it doesn't exist
if not os.path.exists(output_path):
    os.makedirs(output_path)
    print(f"Created directory: {output_path}")

# Get list of .nc files
fluxengine_files = [f for f in os.listdir(fluxengine_data_path) if f.endswith('.nc')]

for file in fluxengine_files:
    full_path = os.path.join(fluxengine_data_path, file)
    
    try:
        # Open dataset
        ds = xr.open_dataset(full_path, decode_times=True)
        print(f"Processing {file}...")
        
        #interpolate onto dan grid
        ds_interp = ds.interp(latitude=lat_grid, longitude=lon_grid, method='linear')
        
        # Save the new file to the specific output folder
        output_filename = "1deg_" + file
        save_path = os.path.join(output_path, output_filename)
        ds_interp.to_netcdf(save_path)

        print(f"Saved to: {save_path}")
        
    except Exception as e:
        print(f"Error processing {file}: {e}")

print("All files processed.")


# Now that the fluxengine output has been downscaled to 1deg on the global monthly
# Dan F data grid I can extract the geographical area covered by the MAXSS data
# from the global monthly data


#Load in a single MAXSS data file to create the mask.
fluxengine_1_deg_1_day_file = 'E:/MAXSS_working_directory/output/MAXSS_RUN/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/1_degree_resolution_output/1deg_MAXSS_2010176N16278_AL012010_ALEX_DOY_160_DATE_2010_06_09.nc'

# Load with xarray (xr) instead of netCDF4 (nc)
fluxengine_1_deg_1_day = xr.open_dataset(fluxengine_1_deg_1_day_file)

#create mask
valid_mask = fluxengine_1_deg_1_day['OF'][0,:,:].notnull()


#Load global data with xarray for easier masking
global_monthly_ds = xr.open_dataset('E:/MAXSS_working_directory/flux_subset_extracted.nc')

# Rename global dims to match the mask dims so they can overlay
global_monthly_ds = global_monthly_ds.rename({'lat': 'latitude', 'lon': 'longitude'})

# Apply the mask to ALL time steps
# Since valid_mask has no 'time' dimension, xarray applies it to every time step automatically.
global_extracted_all_times = global_monthly_ds['flux'].where(valid_mask)

# Save the updated netcdf
# The output file will now contain 3 time steps, but spatially clipped to the storm area.
output_path = "E:/MAXSS_working_directory/global_background_all_3_months.nc"
global_extracted_all_times.to_netcdf(output_path)

print(f"Saved 3-month masked data to: {output_path}")


#sanity check the data
#global monthly average
np.nansum(global_extracted_all_times[:,:,0])

np.nansum(fluxengine_1_deg_1_day['OF'][10,:,:])


## NOW I NEED TO PLOT THE DATA

import xarray as xr
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd

# --- 1. Define the Area Function ---
def get_grid_area_in_m2(latitude_array):
    """Returns grid cell areas (m^2) for a given latitude array."""
    R = 6371000  # Earth radius in meters
    lat_rad = np.radians(latitude_array)
    dy = (np.pi * R) / 180 
    dx = (np.pi * R / 180) * np.cos(lat_rad)
    return dx * dy

# --- 2. Setup Paths ---
fluxengine_data_path = "E:/MAXSS_working_directory/output/MAXSS_RUN/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/"
fluxengine_files = [f for f in os.listdir(fluxengine_data_path) if f.endswith('.nc')]

# Store the full time series here
all_timestamps = []
all_hourly_tg = []  # Changed variable name to Tg

print("Processing files...")

# --- 3. Process Storm Data (Storm Flux) ---
for file in fluxengine_files:
    full_path = os.path.join(fluxengine_data_path, file)
    
    # Open dataset
    ds = xr.open_dataset(full_path)
    
    # Calculate Area
    area_m2 = get_grid_area_in_m2(ds.latitude.values)
    area_da = xr.DataArray(area_m2, coords={'latitude': ds.latitude}, dims='latitude')
    
    # --- CORE CALCULATION ---
    # 1. Multiply Flux * Area -> Total Grams per Day
    flux_rate_total_grams = ds['OF'] * area_da
    
    # 2. Sum over space -> Global Grams per Day
    hourly_rate_global_grams = flux_rate_total_grams.sum(dim=['latitude', 'longitude'])
    
    # 3. Divide by 24 -> Grams per Hour
    actual_mass_hourly_grams = hourly_rate_global_grams / 24.0
    
    # 4. Convert Grams to Teragrams (1e12)
    # 1 Tg = 1,000,000,000,000 grams
    actual_mass_hourly_tg = actual_mass_hourly_grams / 1e12
    
    # Store Data
    all_timestamps.append(ds.time.values)
    all_hourly_tg.append(actual_mass_hourly_tg.values)

# Flatten and prepare Storm Data for plotting
dates_flat = np.concatenate(all_timestamps)
values_flat = np.concatenate(all_hourly_tg)

df_plot = pd.Series(data=values_flat, index=dates_flat)
df_plot = df_plot.sort_index()

# --- 4. Process Global Data (Background Flux) ---
global_monthly_datafile = "E:/MAXSS_working_directory/global_background_all_3_months.nc"
global_monthly_data = xr.open_dataset(global_monthly_datafile)

# Ensure dimensions match (lat -> latitude)
if 'lat' in global_monthly_data.dims:
    global_monthly_data = global_monthly_data.rename({'lat': 'latitude', 'lon': 'longitude'})

bg_dates = []
bg_values = []

# Loop through months
for i in range(len(global_monthly_data.time)):
    date = global_monthly_data.time.values[i]
    global_monthly_co2_flux = global_monthly_data['flux'].isel(time=i)
    
    # Area calculation
    area_m2 = get_grid_area_in_m2(global_monthly_data.latitude.values)
    area_da = xr.DataArray(area_m2, coords={'latitude': global_monthly_data.latitude}, dims='latitude')
    
    # Total Grams per Day
    total_grams_daily = (global_monthly_co2_flux * area_da).sum(dim=['latitude', 'longitude'])
    
    # Convert: (Grams / 24 hours) / 1e12 for Tg
    bg_tg_hourly = (total_grams_daily / 24.0) / 1e12
    
    bg_dates.append(date)
    bg_values.append(bg_tg_hourly.item())

# --- 5. Plotting ---
plt.figure(figsize=(12, 6))

# Plot Storm Flux (Teal)
plt.plot(df_plot.index, df_plot.values, marker='o', markersize=2, 
         linestyle='-', color='teal', label='FluxEngine MAXXS output')

# Plot Background Flux (Orange)
plt.plot(bg_dates, bg_values, color='orange', marker='s', markersize=8, 
         linestyle='--', linewidth=2, label='Global monthly (D.Ford) data')

# Add zero line
plt.axhline(0, color='black', linewidth=0.8)

# Labels
plt.ylabel("Total Carbon Exchange (Tg C hr$^{-1}$)")
plt.xlabel("Date")
plt.title("Hourly Carbon Flux - Storm Alex (2010)")
plt.legend()
plt.grid(True, alpha=0.3)

plt.show()









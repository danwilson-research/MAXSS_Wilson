# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 14:12:31 2026

@author: dw557
"""

import datetime
import netCDF4 as nc
import os

#Location of fluxengine data
fluxengine_data_path = "E:/MAXSS_working_directory/output/MAXSS_RUN/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/"

#Load in the FluxEngine flux data
fluxengine_files = os.listdir(fluxengine_data_path)

## 1: extract the months needed to compare to Dan F data

# Initialize a list to hold the full datetime objects
all_dates = []

# Loop through FluxEngine data to get required dates
for file in fluxengine_files:
    # Open the dataset
    one_day_of_data = nc.Dataset(fluxengine_data_path + file)


    # flux_data = one_day_of_data.variables['OF'][:]
    
    #We assume the time variable exists and the first index [0] is the relevant timestamp
    timestamp_val = float(one_day_of_data.variables['time'][0])
    
    # Convert timestamp to a datetime object
    current_datetime = datetime.datetime.fromtimestamp(timestamp_val, datetime.timezone.utc)

    # Convert to just a date object (removes hours/min/sec)
    just_the_date = current_datetime.date()

    # Add to running list
    all_dates.append(just_the_date)
    
    # Close the file to free up memory
    one_day_of_data.close()

# Check if we found any dates to avoid errors
if all_dates:
    # Get the start (earliest) and end (latest) dates
    start_date = min(all_dates)
    end_date = max(all_dates)

    print(f"Start Date: {start_date}")
    print(f"End Date:   {end_date}")
else:
    print("No dates were found in the file list.")
    

## 2: extract the months of data from Dan F data

import numpy as np

# Location of Dan F global monthly flux data
global_monthly_data_path = 'E:/MAXSS_working_directory/Fordetal_UExP-FNN-U_surface-carbonate-system_v2025-1.nc'
global_monthly_flux = nc.Dataset(global_monthly_data_path)
time_var = global_monthly_flux.variables['time']

# Identify which (Year, Month) pairs we need 
# We use a set to automatically handle duplicates (so we only list 'June 2010' once)
needed_year_months = set((d.year, d.month) for d in all_dates)

print(f"Looking for data matching these months (Year, Month): {sorted(list(needed_year_months))}")

# Convert global time numbers to date objects ---
# nc.num2date automatically handles the "days since 1970-01-15" logic using the file's units
global_dates = nc.num2date(time_var[:], units=time_var.units)

# Find the indices in the global file ---
indices_to_extract = []

for index, dt_obj in enumerate(global_dates):
    # Check if the global file's (Year, Month) matches what we need
    # We access .year and .month directly to avoid cftime errors
    if (dt_obj.year, dt_obj.month) in needed_year_months:
        indices_to_extract.append(index)

# Extract the data ---
if indices_to_extract:
    # We generally want a contiguous slice (e.g., from first match to last match)
    start_index = min(indices_to_extract)
    end_index = max(indices_to_extract) + 1 # +1 for python slicing
    
    print(f"Found matching global data! Slicing from index {start_index} to {end_index}")
    
    #Extract flux data
    global_flux_subset = global_monthly_flux.variables['flux'][:,:,start_index:end_index] 


#save the extracted data into a new much smaller netcdf (code assisatnce from Gemini)

# 1. Define the output filename
output_filename = "flux_subset_extracted.nc"

print(f"Creating new NetCDF file: {output_filename}...")

# 2. Open the new file in 'write' mode
with nc.Dataset(output_filename, "w", format="NETCDF4") as new_nc:
    
    # --- A. SETUP DIMENSIONS ---
    # We infer the sizes from the shape of the data you just extracted.
    # Based on your previous snippet, shape is likely (360, 180, Time) -> (Lon, Lat, Time)
    n_lon = global_flux_subset.shape[0]
    n_lat = global_flux_subset.shape[1]
    n_time = global_flux_subset.shape[2]

    # Create dimensions in the new file
    new_nc.createDimension("lon", n_lon)
    new_nc.createDimension("lat", n_lat)
    new_nc.createDimension("time", n_time) # The subset size

    # --- B. CREATE & FILL COORDINATE VARIABLES ---
    
    # 1. Longitude (Copy full array from source)
    out_lon = new_nc.createVariable("lon", "f4", ("lon",))
    # Try to find the correct name in the source file (usually 'lon' or 'longitude')
    src_lon_name = 'lon' if 'lon' in global_monthly_flux.variables else 'longitude'
    if src_lon_name in global_monthly_flux.variables:
        src_lon = global_monthly_flux.variables[src_lon_name]
        out_lon.units = getattr(src_lon, 'units', 'degrees_east')
        out_lon[:] = src_lon[:] # Copy all longitudes

    # 2. Latitude (Copy full array from source)
    out_lat = new_nc.createVariable("lat", "f4", ("lat",))
    src_lat_name = 'lat' if 'lat' in global_monthly_flux.variables else 'latitude'
    if src_lat_name in global_monthly_flux.variables:
        src_lat = global_monthly_flux.variables[src_lat_name]
        out_lat.units = getattr(src_lat, 'units', 'degrees_north')
        out_lat[:] = src_lat[:] # Copy all latitudes

    # 3. Time (CRITICAL: Copy only the SLICED time values)
    out_time = new_nc.createVariable("time", "f4", ("time",))
    # We already have 'time_var' defined in your script
    out_time.units = getattr(time_var, 'units', 'unknown')
    out_time.calendar = getattr(time_var, 'calendar', 'standard')
    # Use the same start/end index you used for the flux data
    out_time[:] = time_var[start_index:end_index]

    # --- C. SAVE THE FLUX DATA ---
    # We create the variable with dimensions (lon, lat, time) to match your array shape
    out_flux = new_nc.createVariable("flux", "f4", ("lon", "lat", "time"), fill_value=-9999)
    
    # Copy units from original flux variable
    src_flux = global_monthly_flux.variables['flux']
    out_flux.units = getattr(src_flux, 'units', '')
    out_flux.long_name = getattr(src_flux, 'long_name', 'Surface Flux')

    # Write the subset data you extracted in your previous step
    out_flux[:] = global_flux_subset

print("Success! File saved.")

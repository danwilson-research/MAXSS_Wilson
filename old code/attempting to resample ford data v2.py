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



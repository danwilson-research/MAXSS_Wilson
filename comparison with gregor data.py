# -*- coding: utf-8 -*-
"""
Created on Fri Feb  6 10:28:12 2026

@author: dw557
"""

# Script to compare the results from my fluxengine run for Hurricane MARIA, with
# output from the gregor paper OceanSODAETHZv2

#Import required packages
import os
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

#Define target comparison location and data storage location
target_lon = -73 #-72 #
target_lat = 30 #29 # 
fluxengine_output_directory = ("E:/MAXSS_working_directory/output/MAXSS_RUN/maxss/storm-atlas/ibtracs/north-atlantic/2017/2017260N12310_AL152017_MARIA")

#Load in data from my fluxengine output
fluxengine_files = os.listdir(fluxengine_output_directory)

#Section of code to find closest value to target lon and lat in fluxengine output

#load in one file
fluxengine_output = nc.Dataset(fluxengine_output_directory + '/' + fluxengine_files[0])

fluxengine_lons = fluxengine_output['longitude'][:]
fluxengine_lats = fluxengine_output['latitude'][:]

# Find the index of the latitude closest to the target
lat_idx = (np.abs(fluxengine_lats - target_lat)).argmin()

# Find the index of the longitude closest to the target
lon_idx = (np.abs(fluxengine_lons - target_lon)).argmin()

# Sanity check: Print the actual coordinates found to ensure they match closely
found_lat = fluxengine_lats[lat_idx]
found_lon = fluxengine_lons[lon_idx]

print(f"Target Location: {target_lat}, {target_lon}")
print(f"Closest FluxEngine Grid Point: {found_lat:.4f}, {found_lon:.4f}")
print(f"Indices for Fluxengine: Lat index={lat_idx}, Lon index={lon_idx}")

# Initialize lists to store the timeseries data
timeseries_flux = []
timeseries_time = []

# Iterate through files
for file in fluxengine_files:
    # Skip non-NetCDF files (e.g. hidden system files)
    if not file.endswith('.nc'):
        continue
        
    # Load in single file
    file_path = os.path.join(fluxengine_output_directory, file)
    
    #use with loop to ensure automatic closure of netcdf when done with it
    with nc.Dataset(file_path, 'r') as ds:
        # Get the length of the time dimension (should be 24)
        time_len = ds.variables['time'].shape[0]

        # Extract 'OF' at the target location for ALL time steps in this file
        # Slicing [:, lat_idx, lon_idx] gets all times at that specific lat/lon
        flux_vals = ds.variables['OF'][:, lat_idx, lon_idx]
        
        # Append these values to the main list
        # We use extend() because flux_vals is a list/array of 24 items
        timeseries_flux.extend(flux_vals)

        # --- Extract and Convert Time ---
        raw_times = ds.variables['time'][:]
        time_units = ds.variables['time'].units
        
        try:
            time_calendar = ds.variables['time'].calendar
        except AttributeError:
            time_calendar = 'standard'

        # Convert the array of 24 raw time numbers to datetime objects
        dates = nc.num2date(raw_times, units=time_units, calendar=time_calendar)
        
        # Extend the master time list
        timeseries_time.extend(dates)

print(f"Total time points extracted: {len(timeseries_flux)}") 
# Result should be (number of files * 24)
import datetime

# Convert cftime objects to standard Python datetime objects
# We create a new list 'timeseries_time_fixed'
timeseries_time_fixed = [
    datetime.datetime(d.year, d.month, d.day, d.hour, d.minute, d.second) 
    for d in timeseries_time]

# Path to OceanSODA file
oceansoda_path = "E:/MAXSS_working_directory/fgco2-2010s-8D_25km-OceanSODAETHZv2.2024r01.nc" 
oceansoda_var_name = 'fgco2' # variable name for flux

# --- Load OceanSODA Data ---
print(f"Loading OceanSODA data from: {oceansoda_path}")
os_ds = nc.Dataset(oceansoda_path, 'r')

# 1. Extract Coordinate Arrays
os_lats = os_ds.variables['lat'][:]
os_lons = os_ds.variables['lon'][:]
os_times = os_ds.variables['time'][:]

# 2. Find closest Lat/Lon indices (Grid is likely different from FluxEngine)
os_lat_idx = (np.abs(os_lats - target_lat)).argmin()
os_lon_idx = (np.abs(os_lons - target_lon)).argmin()

print(f"OceanSODA Closest Grid: {os_lats[os_lat_idx]:.2f}N, {os_lons[os_lon_idx]:.2f}E")

# 3. Extract Time and Convert
# OceanSODA usually uses "seconds since 1980..." or similar
os_time_units = os_ds.variables['time'].units
try:
    os_calendar = os_ds.variables['time'].calendar
except AttributeError:
    os_calendar = 'standard'

os_dates_raw = nc.num2date(os_times, units=os_time_units, calendar=os_calendar)

# Fix cftime objects if necessary (similar to what we did before)
os_dates_fixed = []
for d in os_dates_raw:
    # Handle both standard datetime and cftime objects
    if isinstance(d, datetime.datetime):
        os_dates_fixed.append(d)
    else:
        # If it's a cftime object, convert to standard datetime
        os_dates_fixed.append(datetime.datetime(d.year, d.month, d.day, d.hour, d.minute))

# 4. Extract Flux Data
# Slicing: [:, os_lat_idx, os_lon_idx] gets all time steps for that location
os_flux = os_ds.variables[oceansoda_var_name][:, os_lat_idx, os_lon_idx]

print(f"Extracted {len(os_flux)} time points from OceanSODA.")

#5. string format for lon/lat
lat_str = f"{abs(target_lat)}°{'N' if target_lat >= 0 else 'S'}"
lon_str = f"{abs(target_lon)}°{'E' if target_lon >= 0 else 'W'}" 

# Convert OceanSODA (mmol) to match FluxEngine (g)
# Factor: 1 mmol C = 0.012011 g C
conversion_factor = 0.012011
os_flux_gC = os_flux * conversion_factor

# --- 2. Plotting ---
plt.figure(figsize=(12, 6))

# Plot FluxEngine (The Storm)
plt.plot(timeseries_time_fixed, timeseries_flux, 
         label='FluxEngine (Hourly Storm Run)', 
         color='blue', linewidth=1,zorder=2)

# Plot OceanSODA (The Background)
# We plot the full series, but we will limit the X-axis view next
plt.step(os_dates_fixed, os_flux_gC, 
         label='OceanSODA-ETHZv2 (8-Day Step)', 
         color='red', linewidth=1, where='mid',zorder=1) # mid used as dataset uses center-adjusted timestamps

# --- 3. Formatting ---
# Set the viewing window to match your FluxEngine run exactly
start_date = timeseries_time_fixed[0]
end_date = timeseries_time_fixed[-1]

# Add a small buffer (e.g., 2 days) to the view so lines don't touch the edges
plt.xlim(start_date - datetime.timedelta(days=2), 
         end_date + datetime.timedelta(days=2))

plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
plt.gcf().autofmt_xdate()

plt.ylabel(r'Air-Sea Flux ($g\ C\ m^{-2}\ day^{-1}$)')
plt.xlabel('Date')
plt.title(f' Hurricane Maria Flux Comparison: ({lat_str}, {lon_str})\nHourly FluxEngine vs 8 Day OceanSODAETHZv2')
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()

plt.show()






# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 12:28:19 2026

@author: dw557
"""

#Install required packages
import os
from os import path
import sys
import gc
import warnings
from glob import glob
import netCDF4 as nc
from string import Template;
import numpy as np;
import calendar;
import pandas as pd
from scipy.interpolate import griddata
from netCDF4 import num2date, Dataset
from scipy.spatial import cKDTree
import itertools
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import xarray as xr # Using xarray makes handling these cubes much easier

MAXSS_working_directory = "E:/MAXSS_working_directory"; 
    
#note to use the same file structure used by the project r.g. #maxss/storm-atlas/ibtracts/region/year/storm
os.chdir(MAXSS_working_directory);
print("Working directory is now:", os.getcwd());

#list of all the basins in MAXSS 
#MAXSS_regions=["east-pacific","north-atlantic","north-indian","south-atlantic","south-indian","south-pacific","west-pacific"]
MAXSS_regions=["north-atlantic"]

#### Loop through the regions in MAXSS storm dataset
for region in MAXSS_regions:
    
    #define the directory for the region
    region_directory = path.join(MAXSS_working_directory+"\\maxss\\storm-atlas\\ibtracs\\{0}".format(region));
    
    #look for the year subfolders

    #get a list of the years
    year_list = []
    for entry_name in os.listdir(region_directory):
        entry_path = path.join(region_directory, entry_name)
        if path.isdir(entry_path):
            year_list.append(entry_name)
    #get a list of the paths for each year folder
    year_directory_list=glob(region_directory+"/*/", recursive = True)
    
    #define to loop through years
    MAXSS_years=year_list
    
#### Loop through the years in the MAXSS storm dataset 
    year_counter=0
    for year in MAXSS_years:
        #get a list of the storms
        storm_list = []
        for entry_name in os.listdir(year_directory_list[year_counter]):
            entry_path = path.join(year_directory_list[year_counter], entry_name)
            if path.isdir(entry_path):
                storm_list.append(entry_name)
        #get a list of the paths for each year folder
        storm_directory_list=glob(year_directory_list[year_counter]+"/*/", recursive = True)
        year_counter=year_counter+1
        
        #define to loop through years
        MAXSS_storms=storm_list
        
        storm_counter=0
        #### Loop through the storms for each year in the MAXSS storm dataset 
        for storm in MAXSS_storms:
            
            ## --- REMOVE SECTION ONCE TESTING COMPLETE --- ##
            if any(name in storm for name in [ ]): #"BONNIE", "COLIN", "MARIA", "RINA" 
                print(f"Skipping storm: {storm}")
                storm_counter += 1 # Important: increment the counter before skipping
                continue

            #directory for storm being processes
            storm_dir=storm_directory_list[storm_counter]
            #directory for storm being processes relative to current working directory
            storm_dir_relative = path.join("\\maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}".format(region,year,storm));
            #output directory 
            output_location=path.join(MAXSS_working_directory+"\\output\\maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}".format(region,year,storm))
            
            #### regrid data for this storm
                          
            # need to get the identifier from the storm name as it is used in .nc file string
            storm_id=storm.split('_')[1]
            
            if region=="north-atlantic":
                region_id="NA"
            elif region=="east-pacific":
                region_id="EP"
            elif region=="north-indian":
                region_id="EP"
            elif region=="south-atlantic":
                region_id="SA"
            elif region=="south-indian":
                region_id="SI"
            elif region=="south-pacific":
                region_id="SP"
            elif region=="west-pacific":
                region_id="WP"  
                
            #### MAXSS L4 Wind data
            #everything is being resampled to wind data resolution
            #we need to extract information from it though
            #we also need to calculate wind speed from east and west components
            
            #### load data
            winds_nc = nc.Dataset(path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\MAXSS_HIST_TC_{3}_{1}_{4}_MAXSS_HIST_TC_L4.nc".format(region,year,storm,region_id,storm_id)));
            
            # 1. Extract the required variables from the NetCDF file
            # We use [:] to load the data from the NetCDF object into a NumPy array
            spatial_mask = winds_nc.variables['data_spatial_mask'][:] 
            track_times = winds_nc.variables['__track_time'][:]
            track_time_units = winds_nc.variables['__track_time'].units # Needed for datetime conversion
            
            # 2. Find out if the storm EVER passed over each grid cell
            # We check if the maximum value along the time axis (axis 0) is 1
            ever_in_storm_mask = np.max(spatial_mask, axis=0) == 1
            
            # 3. Find the index of the FIRST time the storm hit each grid cell
            # argmax() returns the index of the first occurrence of the maximum value (which is 1)
            first_hit_indices = np.argmax(spatial_mask, axis=0)
            
            # 4. Create an empty 2D array filled with NaNs to hold our arrival times
            # It will have the same shape as your lat/lon grid
            grid_shape = (winds_nc.variables['lat'].size, winds_nc.variables['lon'].size)
            storm_arrival_times_numeric = np.full(grid_shape, np.nan)
            
            # 5. Map the actual track times to the grid, but only where the storm actually hit
            storm_arrival_times_numeric[ever_in_storm_mask] = track_times[first_hit_indices[ever_in_storm_mask]]
            
            # 6. (Optional but recommended) Convert the numeric times to standard Python datetimes
            # We flatten the array to convert the valid times, then reshape it back
            storm_arrival_datetimes = np.full(grid_shape, pd.NaT) # Fill with "Not a Time" initially
            
            # Extract only the valid numeric times where the storm hit
            valid_numeric_times = storm_arrival_times_numeric[ever_in_storm_mask]
            
            # Convert them to datetime objects
            valid_datetimes = num2date(valid_numeric_times, units=track_time_units, 
                                       calendar=winds_nc.variables['__track_time'].calendar if hasattr(winds_nc.variables['__track_time'], 'calendar') else 'standard')
            
            # Put the datetimes back into our 2D grid
            storm_arrival_datetimes[ever_in_storm_mask] = valid_datetimes

            print(f"Calculated storm arrival times for {storm}.")
            
            ## Calculate the departure time of the storm in each grid cell
            
            # 1. Get the total number of time steps in our track
            num_time_steps = spatial_mask.shape[0]
            
            # 2. Reverse the spatial mask along the time axis (axis 0)
            # [::-1] is Python's slicing syntax to step backwards
            reversed_mask = spatial_mask[::-1, :, :]
            
            # 3. Find the index of the "first" hit in our reversed time array
            last_hit_indices_reversed = np.argmax(reversed_mask, axis=0)
            
            # 4. Convert the reversed indices back into normal, forward-moving indices
            last_hit_indices = (num_time_steps - 1) - last_hit_indices_reversed
            
            # 5. Create an empty 2D array to hold our departure times
            storm_departure_times_numeric = np.full(grid_shape, np.nan)
            
            # 6. Map the actual track times to the grid for the departure moments
            storm_departure_times_numeric[ever_in_storm_mask] = track_times[last_hit_indices[ever_in_storm_mask]]
            
            # 7. (Optional) Convert to standard Python datetimes, just like we did for arrival
            storm_departure_datetimes = np.full(grid_shape, pd.NaT)
            valid_departure_numeric_times = storm_departure_times_numeric[ever_in_storm_mask]
            
            valid_departure_datetimes = num2date(valid_departure_numeric_times, units=track_time_units, 
                                       calendar=winds_nc.variables['__track_time'].calendar if hasattr(winds_nc.variables['__track_time'], 'calendar') else 'standard')
            
            storm_departure_datetimes[ever_in_storm_mask] = valid_departure_datetimes
            
            print(f"Calculated storm departure times for {storm}.")
            
            ## Calculate if any cells hit twice by storm?
            
            # 1. Calculate the total number of hours the storm was over each cell
            # Summing the 1s and 0s along the time axis (axis 0) gives us the total hours
            total_hours_in_storm = np.sum(spatial_mask, axis=0)
            
            # 2. Calculate the "window" of time from first arrival to last departure
            # We add 1 because if it arrived and left in the same hour, the duration is 1, not 0
            storm_window_duration = (last_hit_indices - first_hit_indices) + 1
            
            # 3. Create a boolean mask highlighting areas hit multiple times
            # It checks if the actual hours are less than the window duration, AND ensures 
            # we are only looking at cells that were actually in the storm to begin with.
            hit_multiple_times_mask = (total_hours_in_storm < storm_window_duration) & ever_in_storm_mask
            
            # 4. Count how many grid cells this happened to and print it for your records
            number_of_cells_hit_multiple_times = np.sum(hit_multiple_times_mask)
            print(f"Grid cells hit multiple times by {storm}: {number_of_cells_hit_multiple_times}")
            
            ## TEST PLOT ##
            
            # 1. Extract latitude and longitude arrays
            lats = winds_nc.variables['lat'][:]
            lons = winds_nc.variables['lon'][:]
            track_lats = winds_nc.variables['__track_lat'][:]
            track_lons = winds_nc.variables['__track_lon'][:]
            
            # 2. Create a meshgrid for plotting
            lon_grid, lat_grid = np.meshgrid(lons, lats)
            
            # 3. Set up the figure and the map projection
            fig = plt.figure(figsize=(10, 8))
            ax = plt.axes(projection=ccrs.PlateCarree())
            
            # 4. Add geographical features
            ax.add_feature(cfeature.COASTLINE, linewidth=1.5)
            ax.add_feature(cfeature.BORDERS, linestyle=':')
            ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
            ax.add_feature(cfeature.OCEAN, facecolor='azure', zorder=0)
            
            # 5a. Plot the BASE LAYER: Arrival Times
            base_mesh = ax.pcolormesh(lon_grid, lat_grid, storm_arrival_times_numeric, 
                                      transform=ccrs.PlateCarree(), cmap='viridis', shading='auto', zorder=1)
            
            # 5b. Plot the OVERLAY LAYER: Multiple Hits (in Red)
            red_overlay_data = np.where(hit_multiple_times_mask, 1.0, np.nan)
            red_cmap = mcolors.ListedColormap(['red'])
            
            ax.pcolormesh(lon_grid, lat_grid, red_overlay_data, 
                          transform=ccrs.PlateCarree(), cmap=red_cmap, shading='auto', zorder=2)
            
            # 5c. Plot the STORM TRACK
            # Notice we assign the plot to a variable 'track_line' so we can pass it to the legend
            track_line, = ax.plot(track_lons, track_lats, color='orange', linewidth=2, marker='o', markersize=3, 
                                  transform=ccrs.PlateCarree(), zorder=3, label='Storm Center Track')
            
            # ---> NEW: Create a custom patch for the legend
            red_patch = mpatches.Patch(color='red', label='Multiple Hits (Gaps)')
            
            # ---> NEW: Tell the legend to use both the track line and our custom red patch
            ax.legend(handles=[track_line, red_patch], loc='lower right')
            
            # 6. Add a colorbar linked ONLY to the base arrival time layer
            cbar = plt.colorbar(base_mesh, ax=ax, orientation='vertical', shrink=0.7)
            cbar.set_label(f"Arrival Time ({track_time_units})")
            
            # Focus the map exactly on the bounding box
            ax.set_extent([np.min(lons), np.max(lons), np.min(lats), np.max(lats)], crs=ccrs.PlateCarree())
            
            plt.title(f"Arrival Times & Multiple Hits for {storm}")
            
            plt.show()
            
    

            # Remember to close your NetCDF file at the end of the loop to free up memory!
            winds_nc.close()
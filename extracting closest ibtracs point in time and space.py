# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 12:28:19 2026

@author: dw557
"""

#Initially this script loads in the ibtracs data and extracts the data for each releveant storm

#Load required packages
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import pandas as pd
from os import path
from glob import glob

# change to plt.close from plt.show when running for many storms

# --- CONFIGURATION ---
MAXSS_working_directory = "E:/MAXSS_working_directory"

#Load in ibtracs dataset
ibtracs_dataset_path = "E:/MAXSS_working_directory/ibtracs.ALL.list.v04r01.csv"
ibtracs_data = pd.read_csv(ibtracs_dataset_path, skiprows=[1], low_memory=False)

#CHANGE WHEN I MOVE TO INLCUDING MORE REGIONS
MAXSS_regions = ["north-atlantic"]

# --- MAIN LOOP ---
for region in MAXSS_regions:
    
    # Define the directory for the region
    region_directory = path.join(MAXSS_working_directory, "output/MAXSS_RUN/maxss/storm-atlas/ibtracs", region)
    
    # Get all year directories
    year_directory_list = glob(path.join(region_directory, "*/"))
    
    for year_dir in year_directory_list:
        # Extract year name from path for labelling if needed
        year_name = os.path.basename(os.path.normpath(year_dir))
        
        # Get list of storms in this year directory
        storm_directory_list = glob(path.join(year_dir, "*/"))

        for storm_path in storm_directory_list:
            #print(storm_path)
        
            # Get the storm name from the folder path
            storm_name = os.path.basename(os.path.normpath(storm_path))
            
            storm_id = storm_name.split('_')[0]
            
            ## --- TESTING: Skip logic --- ##
            # List the storms you DO NOT want to process
            storms_to_skip = []
            #storms_to_skip = ["ALEX", "BONNIE", "MARIA", "RINA"]

            # If the current storm_name DOES contain ANY of the target names, skip it
            if any(name in storm_name for name in storms_to_skip):
                print(f"Skipping storm: {storm_name}")
                continue
            ## --------------------------- ##
            
            print(f'Working on storm: {storm_name}')
            
            # --- 1. Construct the Path to the NetCDF mask ---
            # Based on your example: E:\MAXSS_working_directory\maxss\storm-atlas\ibtracs\north-atlantic\2010\STORM_NAME\filename.nc
            mask_filename = "Resampled_for_fluxengine_MAXSS_land_fraction.nc"
            
            # Use path.join to keep it clean and cross-platform
            mask_path = path.join(MAXSS_working_directory, "maxss", "storm-atlas", "ibtracs", 
                                  region, year_name, storm_name, mask_filename)

            # --- 2. Load the Mask Data ---
            if path.exists(mask_path):
                try:
                    ds_mask = xr.open_dataset(mask_path)
                    
                    # Access the land fraction variable
                    land_frac = ds_mask['land proportion'] 
                    
                    # --- NEW: Create the mask for values between 0 and 1 ---
                    # This creates a DataArray of True/False
                    valid_mask = (land_frac >= 0) & (land_frac <= 1)
                    
                    #Convert to 1s and 0s if you prefer (True=1, False=0)
                    valid_mask = valid_mask.astype(int)

                    print(f"Successfully created mask for {storm_name}")

                    
                except Exception as e:
                    print(f"Error loading mask for {storm_name}: {e}")
            
            # 1. Extract and Clean Data
            storm_data = ibtracs_data[ibtracs_data['SID'] == storm_id].copy()
            
            # CRITICAL: Convert LAT/LON to numeric (they often load as strings)
            storm_data['LAT'] = pd.to_numeric(storm_data['LAT'], errors='coerce')
            storm_data['LON'] = pd.to_numeric(storm_data['LON'], errors='coerce')
            
            # Remove rows with missing coordinates
            storm_data = storm_data.dropna(subset=['LAT', 'LON'])
            
            # 2. Handle Time and Sorting
            if 'ISO_TIME' in storm_data.columns:
                storm_data['ISO_TIME'] = pd.to_datetime(storm_data['ISO_TIME'])
                storm_data = storm_data.sort_values('ISO_TIME')
            
            # 3. Quick Plotting Check
            if not storm_data.empty:
                # 1. Create figure and axis with projection
                fig = plt.figure(figsize=(12, 8))
                ax = plt.axes(projection=ccrs.PlateCarree())
                
                # 2. OVERLAY THE MASK (if it exists)
                if 'valid_mask' in locals() or 'valid_mask' in globals():
                    # We use .where(valid_mask) to make the False values transparent (NaN)
                    # This allows the map features/land to show through the invalid areas
                    mask_to_plot = valid_mask.where(valid_mask == True)
                    
                    mask_to_plot.plot(ax=ax, transform=ccrs.PlateCarree(), 
                                      add_colorbar=False, cmap='Blues', alpha=0.3)

                # 3. Add Map Features
                ax.add_feature(cfeature.LAND, facecolor='lightgray', alpha=0.5)
                ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
                ax.add_feature(cfeature.BORDERS, linestyle=':', alpha=0.4)

                # 4. Plot the IBTrACS Track
                # ax.plot(storm_data['LON'], storm_data['LAT'], 
                #         linestyle='-', color='red', linewidth=1.5, 
                #         label='IBTrACS Track', transform=ccrs.PlateCarree(), zorder=4)
                
                # 5. Plot the IBTrACS Points
                # We use zorder=5 to ensure points are on the very top
                scatter = ax.scatter(storm_data['LON'], storm_data['LAT'], 
                                     c='black',  
                                     s=30, edgecolors='black', linewidth=0.5,
                                     zorder=5, transform=ccrs.PlateCarree(),
                                     label='IBTrACS Points')
            
                # 6. Formatting
                plt.title(f"Storm Track Overlay: {storm_name}\n(SID: {storm_id})")
                ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, alpha=0.3)
                
                plt.legend(loc='lower left')
                
                plt.show()
                
            else:
                print(f"No valid coordinate data found for {storm_id}")
            
            ##NOW EXTRACT THE TIME OF THE CLOSEST APPROACH
            
            # --- Calculate Time of Closest Approach ---

            # --- Calculate Time of Closest Approach ---

            # 1. Extract coordinates as raw numpy arrays for broadcasting
            # Note: Ensure these match the order of dimensions in your mask (usually lat then lon)
            lats = ds_mask.latitude.values  # Shape: (N_lat,)
            lons = ds_mask.longitude.values  # Shape: (N_lon,)
            storm_lats = storm_data['LAT'].values  # Shape: (N_storm,)
            storm_lons = storm_data['LON'].values  # Shape: (N_storm,)
            storm_times = storm_data['ISO_TIME'].values
            
            # 2. Broadcasting for distances
            # We create a (N_storm, N_lat, N_lon) array
            # [:, None, None] expands storm points, [None, :, None] expands lats, [None, None, :] expands lons
            delta_lat = storm_lats[:, np.newaxis, np.newaxis] - lats[np.newaxis, :, np.newaxis]
            delta_lon = storm_lons[:, np.newaxis, np.newaxis] - lons[np.newaxis, np.newaxis, :]
            
            # Squared Euclidean distance (broadcasting handles the expansion)
            dist_sq = delta_lat**2 + delta_lon**2
            
            # 3. Find index of closest storm point for every grid cell
            # This reduces the 3D array back to a 2D array of shape (N_lat, N_lon)
            closest_idx = dist_sq.argmin(axis=0)
            
            # 4. Map indices to times and wrap back into Xarray
            # We use the coordinates and dimensions from the original land_frac/mask
            time_of_closest_approach = xr.DataArray(
                storm_times[closest_idx],
                coords=[ds_mask.latitude, ds_mask.longitude],
                dims=['latitude', 'longitude'],
                name="time_of_closest_approach"
            )
            
            # 5. Apply your valid_mask
            time_of_closest_approach = time_of_closest_approach.where(valid_mask == 1)
            
            # 1. Create the full grid of times BEFORE masking
            time_of_closest_approach_full = xr.DataArray(
                storm_times[closest_idx],
                coords=[ds_mask.latitude, ds_mask.longitude],
                dims=['latitude', 'longitude']
            )
            
            # 2. Calculate the time difference
            start_time = storm_times.min()
            time_delta = time_of_closest_approach_full - start_time
            
            # 3. Convert to float days using the .dt accessor
            # total_seconds() safely converts timedeltas to floats, then we divide by seconds in a day
            time_days = time_delta.dt.total_seconds() / 86400
            
            # 4. NOW apply the mask (the masked areas will become normal float NaNs)
            time_days = time_days.where(valid_mask == 1)
            
            # 5. Plotting
            fig, ax = plt.subplots(figsize=(12, 7), subplot_kw={'projection': ccrs.PlateCarree()})
            
            # The colorbar should now automatically scale from 0 to ~15 days
            mesh = time_days.plot(
                ax=ax, 
                cmap='viridis', 
                add_colorbar=True, 
                cbar_kwargs={'label': 'Days since storm start'}
            )
            
            ax.coastlines()
            plt.title(f"Time of Closest Approach (Relative Days): {storm_name}")
            plt.show()
            
          ### I need to make sure code is neatened and tidied up e.g. comments
          ### I also need to colve the issue of storms going back on themselves
          ### do this with when they first came within x distance buffer
            
            print(f"Calculated closest approach for {storm_name}")
            
            #Memory management
            if 'ds_mask' in locals():
                ds_mask.close()
            
            
            
            
            
            
            
            
            
            
            
            # DYNAMIC PATH ASSIGNMENT (changes with each storm)
            fluxengine_dir = storm_path
            
            # Get list of all NetCDF files in this storm folder
            fluxengine_files = [f for f in os.listdir(fluxengine_dir) if f.endswith('.nc')]
            
            
            
            
            
            if not fluxengine_files:
                print(f"Skipping {storm_name}: No NetCDF files found.")
                continue
            
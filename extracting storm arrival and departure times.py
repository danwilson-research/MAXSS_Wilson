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
import numpy as np;
import pandas as pd
from netCDF4 import num2date, Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import xarray as xr # Using xarray makes handling these cubes much easier
import datetime

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
                
            #### load L4 Wind data
            winds_nc = nc.Dataset(path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\MAXSS_HIST_TC_{3}_{1}_{4}_MAXSS_HIST_TC_L4.nc".format(region,year,storm,region_id,storm_id)));
            
            # 1. Extract the required variables from the NetCDF file
            spatial_mask = winds_nc.variables['data_spatial_mask'][:] 
            track_times = winds_nc.variables['__track_time'][:]
            track_time_units = winds_nc.variables['__track_time'].units # Needed for datetime conversion
            
            # 2. Find out if the storm EVER passed over each grid cell
            # We check if the maximum value along the time axis (axis 0) is 1
            ever_in_storm_mask = np.max(spatial_mask, axis=0) == 1
            
            ## Calculate storm initial arrival time ##
            
            # 1. Find the index of the FIRST time the storm hit each grid cell
            # argmax() returns the index of the first occurrence of the maximum value (which is 1)
            first_hit_indices = np.argmax(spatial_mask, axis=0)
            
            # 2. Create an empty 2D array filled with NaNs to hold our arrival times
            # It has the same shape as lat/lon grid
            grid_shape = (winds_nc.variables['lat'].size, winds_nc.variables['lon'].size)
            storm_arrival_times_numeric = np.full(grid_shape, np.nan)
            
            # 3. Map the actual track times to the grid, but only where the storm actually hit
            storm_arrival_times_numeric[ever_in_storm_mask] = track_times[first_hit_indices[ever_in_storm_mask]]
            
            # 4. Extract only the valid numeric times where the storm hit
            valid_numeric_times = storm_arrival_times_numeric[ever_in_storm_mask]
            
            ## Calculate the departure time of the storm in each grid cell ##
            
            # 1. Get the total number of time steps in track
            num_time_steps = spatial_mask.shape[0]
            
            # 2. Reverse the spatial mask along the time axis (axis 0)
            # [::-1] is Python's slicing syntax to step backwards
            reversed_mask = spatial_mask[::-1, :, :]
            
            # 3. Find the index of the "first" hit in reversed time array
            last_hit_indices_reversed = np.argmax(reversed_mask, axis=0)
            
            # 4. Convert the reversed indices back into normal, forward-moving indices
            last_hit_indices = (num_time_steps - 1) - last_hit_indices_reversed
            
            # 5. Create an empty 2D array to hold departure times
            storm_departure_times_numeric = np.full(grid_shape, np.nan)
            
            # 6. Map the actual track times to the grid for the departure moments
            storm_departure_times_numeric[ever_in_storm_mask] = track_times[last_hit_indices[ever_in_storm_mask]]
            
            #print(f"Calculated storm arrival and departure times for {storm}.")
            
            ## Calculate the start of the pre-storm reference period ##
            
            # 1. Convert numeric arrival times to datetime objects 
            # (These become the special cftime objects to match the NetCDF calendar)
            valid_arrival_datetimes = num2date(valid_numeric_times, units=track_time_units, 
                                               calendar=winds_nc.variables['__track_time'].calendar if hasattr(winds_nc.variables['__track_time'], 'calendar') else 'standard')
            
            # 2. Subtract exactly 15 days using Python's standard timedelta directly!
            # No Pandas conversion needed; cftime objects understand this math perfectly.
            valid_pre_storm_start_times = valid_arrival_datetimes - datetime.timedelta(days=15)
            
            # 3. Create an empty 2D grid to hold our results
            # We use dtype=object and fill it with None so it can safely hold the cftime objects
            pre_storm_ref_period_start = np.full(grid_shape, None, dtype=object)
            
            # 4. Map our calculated 15-day-prior times back onto the exact grid cells the storm hit
            pre_storm_ref_period_start[ever_in_storm_mask] = valid_pre_storm_start_times
            
            # Convert Pre-Storm Start to 'seconds since 1981-01-01'
            # 1. Define exactly what units we want the output to be in
            target_time_units = track_time_units
            
            # 2. Grab the calendar type again to ensure our conversion is perfectly accurate
            storm_calendar = winds_nc.variables['__track_time'].calendar if hasattr(winds_nc.variables['__track_time'], 'calendar') else 'standard'
            
            # 3. Convert the valid cftime objects back into flat numbers
            # We use the nc.date2num function (the exact opposite of num2date)
            valid_pre_storm_start_numeric = nc.date2num(valid_pre_storm_start_times, 
                                                        units=target_time_units, 
                                                        calendar=storm_calendar)
            
            # 4. Create a fresh 2D grid filled with standard NaNs to hold the numbers
            pre_storm_ref_period_start_numeric = np.full(grid_shape, np.nan)
            
            # 5. Map our newly calculated seconds back onto the exact grid cells the storm hit
            pre_storm_ref_period_start_numeric[ever_in_storm_mask] = valid_pre_storm_start_numeric
            
            ## Calculate the end of the Pre-Storm Reference Period (which is 2 days before storm arrival) ##
            
            # 1. Subtract exactly 2 days from the arrival times using Python's standard timedelta
            valid_pre_storm_end_times = valid_arrival_datetimes - datetime.timedelta(days=2)
            
            # 2. Create an empty 2D grid to hold our cftime objects
            pre_storm_ref_period_end = np.full(grid_shape, None, dtype=object)
            
            # 3. Map calculated 2-day-prior times back onto the exact grid cells the storm hit
            pre_storm_ref_period_end[ever_in_storm_mask] = valid_pre_storm_end_times
            
            # 4. Convert the valid cftime objects back into flat numbers
            # We reuse the target_time_units and storm_calendar we defined in the previous step
            valid_pre_storm_end_numeric = nc.date2num(valid_pre_storm_end_times, 
                                                      units=target_time_units, 
                                                      calendar=storm_calendar)
            
            # 5. Create a fresh 2D grid filled with standard NaNs to hold the numbers
            pre_storm_ref_period_end_numeric = np.full(grid_shape, np.nan)
            
            # 6. Map newly calculated seconds back onto the exact grid cells the storm hit
            pre_storm_ref_period_end_numeric[ever_in_storm_mask] = valid_pre_storm_end_numeric
            
            #print(f"Calculated start and end of pre-storm reference period for {storm}.")
            
            ## Calculate 40-day Post-Storm Analysis End ##
            
            # 1. Extract only the valid numeric departure times where the storm hit
            valid_departure_numeric = storm_departure_times_numeric[ever_in_storm_mask]
            
            # 2. Convert numeric departure times to datetime (cftime) objects
            # We reuse the track_time_units and storm_calendar we defined earlier
            valid_departure_datetimes = num2date(valid_departure_numeric, units=track_time_units, 
                                                 calendar=storm_calendar)
            
            # 4. Add exactly 40 days using Python's standard timedelta
            valid_post_storm_end_times = valid_departure_datetimes + datetime.timedelta(days=40)
            
            # 5. Create empty 2D grid to hold cftime objects
            post_storm_analysis_end = np.full(grid_shape, None, dtype=object)
            
            # 6. Map calculated 40-day-post times back onto the exact grid cells the storm hit
            post_storm_analysis_end[ever_in_storm_mask] = valid_post_storm_end_times
            
            # 7. Convert the valid cftime objects back into flat numbers
            valid_post_storm_end_numeric = nc.date2num(valid_post_storm_end_times, 
                                                       units=target_time_units, 
                                                       calendar=storm_calendar)
            
            # 8. Create a fresh 2D grid filled with standard NaNs to hold the numbers
            post_storm_analysis_end_numeric = np.full(grid_shape, np.nan)
            
            # 9. Map our newly calculated seconds back onto the exact grid cells the storm hit
            post_storm_analysis_end_numeric[ever_in_storm_mask] = valid_post_storm_end_numeric
            
            print(f"Calculated storm timings {storm}.")
            
            ## Calculate if any cells hit twice by storm
            
            # 1. Calculate the total number of hours the storm was over each cell
            # Summing the 1s and 0s along the time axis (axis 0) gives the total hours
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
            
            ## --- Create 3D Analysis Period Mask --- ##
            
            # 1. Grab the full data time array (475 steps)
            full_data_times = winds_nc.variables['time'][:]
            full_data_units = winds_nc.variables['time'].units
            
            # 2. Reshape for 3D broadcasting: (475, 1, 1)
            full_data_times_3d = full_data_times[:, np.newaxis, np.newaxis]
            
            # 3. Re-create the 3D mask using the FULL timeline
            # This compares every hour of the data to your spatial 2D arrival/departure maps
            analysis_period_mask_3d = (full_data_times_3d >= pre_storm_ref_period_start_numeric) & \
                                      (full_data_times_3d <= post_storm_analysis_end_numeric)
            
            print(f"New Mask Shape: {analysis_period_mask_3d.shape}") # Should be (475, Lat, Lon)
            
            #Save the storm timings to a netcdf file
            
            #Load in required bits of wind file
            wind_eastward = winds_nc.variables['__eo_eastward_wind'][:]
            wind_northward = winds_nc.variables['__eo_northward_wind'][:]
            
            wind_eastward = wind_eastward.astype('float64')
            wind_northward = wind_northward.astype('float64')
            
            #### get wind data dimensions
            wind_time_dimension=len(wind_eastward)
            wind_lat_dimension=len(wind_eastward[0])
            wind_lon_dimension=len(wind_eastward[0][0])
            
            #### define lat,long,time variables - these are used to build new matrixes
            wind_lat = winds_nc.variables['lat'][:]
            wind_lon = winds_nc.variables['lon'][:]
            wind_time = winds_nc.variables['time'][:]
            wind_dates = num2date(wind_time, winds_nc.variables['time'].units)
             
            #### min and max of wind grid lats and lons   
            min_lat=wind_lat[0]
            max_lat=wind_lat[0-1]                    
            min_lon=wind_lon[0]
            max_lon=wind_lon[0-1] 
            
            # --- Save storm timings AND 3D Mask to a netcdf file ---
            processedFilePath = (path.join(f"maxss\\storm-atlas\\ibtracs\\{region}\\{year}\\{storm}\\Resampled_for_fluxengine_storm_timings_with_mask.nc"))
            ncout = Dataset(processedFilePath, 'w')
            
            # 1. Define Dimensions
            # We now include 'time' which matches your 475-step data timeline
            ncout.createDimension("time", len(full_data_times))
            ncout.createDimension("lat", wind_lat_dimension)
            ncout.createDimension("lon", wind_lon_dimension)
            
            # 2. Create Dimension Variables
            var_time = ncout.createVariable("time", float, ("time",))
            var_time.units = full_data_units
            var_time.calendar = storm_calendar
            var_time[:] = full_data_times

            var_lat = ncout.createVariable("lat", float, ("lat",))
            var_lat.units = "degrees_north"
            var_lat[:] = wind_lat
            
            var_lon = ncout.createVariable("lon", float, ("lon",))
            var_lon.units = "degrees_east"
            var_lon[:] = wind_lon

            # 3. Create the 3D Mask Variable
            # We use 'int8' (or 'byte') to save space since it's just 0s and 1s
            var_mask = ncout.createVariable("analysis_mask", "i1", ("time", "lat", "lon"), 
                                            zlib=True, complevel=4, shuffle=True)
            var_mask.long_name = "3D Analysis Period Mask (-15 to +40 days)"
            var_mask.description = "1 = Valid window for that pixel, 0 = Outside window"
            # We convert the boolean mask to integers (0 and 1) for saving
            var_mask[:] = analysis_period_mask_3d.astype(np.int8)

            # 4. Create the 2D Timing Variables (Same as your original code)
            var_pre_start = ncout.createVariable("pre_storm_start", float, ("lat", "lon"), zlib=True)
            var_pre_start.units = track_time_units
            var_pre_start[:] = pre_storm_ref_period_start_numeric
            
            var_pre_end = ncout.createVariable("pre_storm_end", float, ("lat", "lon"), zlib=True)
            var_pre_end.units = track_time_units
            var_pre_end[:] = pre_storm_ref_period_end_numeric
            
            var_arrive = ncout.createVariable("storm_arrival", float, ("lat", "lon"), zlib=True)
            var_arrive.units = track_time_units
            var_arrive[:] = storm_arrival_times_numeric
            
            var_depart = ncout.createVariable("storm_departure", float, ("lat", "lon"), zlib=True)
            var_depart.units = track_time_units
            var_depart[:] = storm_departure_times_numeric
            
            var_post_end = ncout.createVariable("post_storm_end", float, ("lat", "lon"), zlib=True)
            var_post_end.units = track_time_units
            var_post_end[:] = post_storm_analysis_end_numeric

            ncout.close()
            print(f"Successfully saved 3D mask and timing variables for {storm}.")
            
            #####################################################################################
            # ## --- Sanity Check: Plot the 3D Mask at 4 Timesteps --- ##
            
            # # 1. Extract lat/lon/time arrays
            # lats = winds_nc.variables['lat'][:]
            # lons = winds_nc.variables['lon'][:]
            # # Note: Using the master data time array here
            # data_times = winds_nc.variables['time'][:] 
            
            # # 2. Calculate indices spread evenly through the ENTIRE file
            # # This allows us to see the mask state at the very beginning and very end of the data
            # total_steps = len(full_data_times)
            # time_indices = [0, int(total_steps*0.33), int(total_steps*0.66), total_steps-1]
            
            # fig, axes = plt.subplots(2, 2, figsize=(14, 10), 
            #                          subplot_kw={'projection': ccrs.PlateCarree()})
            # axes = axes.flatten()
            
            # for i, t_idx in enumerate(time_indices):
            #     ax = axes[i]
                
            #     # This will no longer throw an IndexError because axis 0 is now 475
            #     mask_slice = analysis_period_mask_3d[t_idx, :, :].astype(int)
                
            #     mesh = ax.pcolormesh(lons, lats, mask_slice, cmap='coolwarm', 
            #                          vmin=0, vmax=1, transform=ccrs.PlateCarree())
                
            #     ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
                
            #     # Formatting title with correct date
            #     dt_obj = num2date(full_data_times[t_idx], units=full_data_units, calendar=storm_calendar)
            #     ax.set_title(f"Data Step {t_idx}: {dt_obj.strftime('%Y-%m-%d %H:%M')}")
            
            # plt.show()
            
            # import matplotlib.animation as animation

            # # --- 1. Setup the Figure ---
            # fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={'projection': ccrs.PlateCarree()})
            # ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
            # ax.add_feature(cfeature.BORDERS, linestyle=':', alpha=0.5)
            
            # # Initialize the plot with the first frame (Step 0)
            # mask_slice = analysis_period_mask_3d[0, :, :].astype(int)
            # mesh = ax.pcolormesh(lons, lats, mask_slice, cmap='coolwarm', vmin=0, vmax=1, transform=ccrs.PlateCarree())
            
            # # Setup the title and colorbar
            # title = ax.set_title(f"Mask Progression: {storm} | Step 0")
            # cbar = fig.colorbar(mesh, orientation='horizontal', shrink=0.5, pad=0.05)
            # cbar.set_ticks([0.25, 0.75])
            # cbar.set_ticklabels(['Masked', 'Valid Window'])
            
            # # --- 2. Define the Update Function ---
            # def update(frame):
            #     # Update the data in the existing mesh object
            #     mask_slice = analysis_period_mask_3d[frame, :, :].astype(int)
            #     mesh.set_array(mask_slice.flatten())
                
            #     # Update the title with the current date
            #     dt_obj = num2date(full_data_times[frame], units=full_data_units, calendar=storm_calendar)
            #     title.set_text(f"Mask Progression: {storm}\nDate: {dt_obj.strftime('%Y-%m-%d %H:%M')}")
            #     return mesh, title
            
            # # --- 3. Create the Animation ---
            # # We use a step to make the file size manageable; set to 1 for maximum smoothness
            # frame_step = 10 
            # frames_to_run = range(0, len(full_data_times), frame_step)
            
            # ani = animation.FuncAnimation(fig, update, frames=frames_to_run, interval=50, blit=True)
                      
            # # To show in your environment:
            # plt.show()
            
            # # Save the animation as a GIF
            # ani.save(f'{storm}_mask_animation.gif', writer='pillow', fps=15)
            # print(f"Animation saved as {storm}_mask_animation.gif in your working directory.")

################################################################################
           
            ## Test plot showing storm timings and 
            
            # 1. Find the grid cell closest to the true middle of the storm track
            track_lats = winds_nc.variables['__track_lat'][:]
            track_lons = winds_nc.variables['__track_lon'][:]
            grid_lats = winds_nc.variables['lat'][:]
            grid_lons = winds_nc.variables['lon'][:]
            
            # 1. Strip out any NaNs so we are only looking at real track coordinates
            valid_track_mask = ~np.isnan(track_lats) & ~np.isnan(track_lons)
            valid_lats = track_lats[valid_track_mask]
            valid_lons = track_lons[valid_track_mask]
            
            if len(valid_lats) > 0:
                # Get the geographic coordinate right in the middle of the valid track
                mid_track_idx = len(valid_lats) // 2
                center_lat = valid_lats[mid_track_idx]
                center_lon = valid_lons[mid_track_idx]
                
                # 2. The Longitude Fix: Force both track and grid to a 0-360 scale 
                center_lon_360 = center_lon % 360
                grid_lons_360 = grid_lons[:] % 360
                
                # Find the closest matching index in our NetCDF spatial grid
                test_lat_idx = np.abs(grid_lats - center_lat).argmin()
                test_lon_idx = np.abs(grid_lons_360 - center_lon_360).argmin()
                
                # --- NEW: PRINT THE SANITY CHECK ---
                matched_lat = grid_lats[test_lat_idx]
                matched_lon = grid_lons[test_lon_idx]
                #print(f"Aiming for Track Center: Lat {center_lat:.2f}, Lon {center_lon:.2f}")
                #print(f"Snapped to Grid Cell: Lat {matched_lat:.2f}, Lon {matched_lon:.2f} (Indices: {test_lat_idx}, {test_lon_idx})")
                # -----------------------------------
                
                # 3. The Mask Fix: Check if this ideal center cell actually has valid storm data!
                if not ever_in_storm_mask[test_lat_idx, test_lon_idx]:
                    print("Warning: Center track point missed the valid mask. Falling back to a known valid cell...")
                    valid_lat_indices, valid_lon_indices = np.where(ever_in_storm_mask)
                    mid_idx = len(valid_lat_indices) // 2
                    test_lat = valid_lat_indices[mid_idx]
                    test_lon = valid_lon_indices[mid_idx]
                else:
                    # If the mask check passed, finalize our test indices
                    test_lat = test_lat_idx
                    test_lon = test_lon_idx
        
                # 2. Extract the time series for this specific cell
                # Grab the U (eastward) and V (northward) components
                u_ts = winds_nc.variables['__eo_eastward_wind'][:, test_lat, test_lon]
                v_ts = winds_nc.variables['__eo_northward_wind'][:, test_lat, test_lon]
                
                # Calculate the total wind speed magnitude
                wind_ts = np.sqrt(u_ts**2 + v_ts**2)
                
                # Calculate the total wind speed magnitude
                wind_ts = np.sqrt(u_ts**2 + v_ts**2)
        
                # 3. Get the calculated numeric times for this specific cell
                t_pre_start_num = pre_storm_ref_period_start_numeric[test_lat, test_lon]
                t_pre_end_num = pre_storm_ref_period_end_numeric[test_lat, test_lon]
                t_arrive_num = storm_arrival_times_numeric[test_lat, test_lon]
                t_depart_num = storm_departure_times_numeric[test_lat, test_lon]
                t_post_end_num = post_storm_analysis_end_numeric[test_lat, test_lon]
        
                # 4. We need the FULL time axis for the wind, not just the track times!
                full_times = winds_nc.variables['time'][:]
                full_time_units = winds_nc.variables['time'].units
                full_calendar = winds_nc.variables['time'].calendar if hasattr(winds_nc.variables['time'], 'calendar') else 'standard'
                
                # Convert numbers to cftime objects
                raw_plot_times = num2date(full_times, units=full_time_units, calendar=full_calendar)
                
                # FIX: Convert cftime objects into standard Python datetime objects for Matplotlib!
                plot_times = [datetime.datetime(d.year, d.month, d.day, d.hour, d.minute, d.second) for d in raw_plot_times]
                
                # Helper function to safely convert boundary points back to standard datetimes for plotting
                def safe_num2date(num_val):
                    if not np.isnan(num_val):
                        d = num2date(num_val, units=track_time_units, calendar=storm_calendar)
                        # FIX: Convert to standard datetime
                        return datetime.datetime(d.year, d.month, d.day, d.hour, d.minute, d.second)
                    return None
        
                t_pre_start_dt = safe_num2date(t_pre_start_num)
                t_pre_end_dt = safe_num2date(t_pre_end_num)
                t_arrive_dt = safe_num2date(t_arrive_num)
                t_depart_dt = safe_num2date(t_depart_num)
                t_post_end_dt = safe_num2date(t_post_end_num)
        
                # 5. Create the plot!
                plt.figure(figsize=(14, 6))
                
                # Plot the actual wind speed against the FULL timeline
                plt.plot(plot_times, wind_ts, label='Wind Speed', color='black', linewidth=1.5)
                
                # Shade the background using the arrival and departure datetimes directly
                if t_arrive_dt and t_depart_dt:
                    plt.axvspan(t_arrive_dt, t_depart_dt, color='grey', alpha=0.3, label='Active Storm Window')
        
                # Plot the timeline boundaries
                if t_pre_start_dt: plt.axvline(t_pre_start_dt, color='blue', linestyle='--', linewidth=2, label='Pre-Storm Start (-15d)')
                if t_pre_end_dt: plt.axvline(t_pre_end_dt, color='cyan', linestyle='--', linewidth=2, label='Pre-Storm End (-2d)')
                if t_arrive_dt: plt.axvline(t_arrive_dt, color='red', linestyle='-', linewidth=2.5, label='First Arrival')
                if t_depart_dt: plt.axvline(t_depart_dt, color='purple', linestyle='-', linewidth=2.5, label='Final Departure')
                if t_post_end_dt: plt.axvline(t_post_end_dt, color='green', linestyle='--', linewidth=2, label='Post-Storm End (+40d)')
        
                # Formatting
                plt.title(f"Diagnostic Timeline for {storm} (Lat: {matched_lat}, Lon: {matched_lon})")
                plt.ylabel("Wind Speed (m/s)")
                plt.xlabel("Date")
                
                # Move legend outside the plot so it doesn't block the data
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                plt.grid(True, alpha=0.3)
                plt.tight_layout()
                plt.show()
                
            else:
                print(f"No valid storm cells found to plot for {storm}.")
                
            ## Plot of storm track and storm analysis area ##
            
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
            
            # 6. Create a custom patch for the legend
            red_patch = mpatches.Patch(color='red', label='Multiple Hits (Gaps)')
            
            # 7. Tell the legend to use both the track line and our custom red patch
            ax.legend(handles=[track_line, red_patch], loc='lower right')
            
            # 8. Add a colorbar linked ONLY to the base arrival time layer
            cbar = plt.colorbar(base_mesh, ax=ax, orientation='vertical', shrink=0.7)
            cbar.set_label(f"Arrival Time ({track_time_units})")
            
            # Focus the map exactly on the bounding box
            ax.set_extent([np.min(lons), np.max(lons), np.min(lats), np.max(lats)], crs=ccrs.PlateCarree())
            
            # draw_labels=True gives you the Lat/Lon coordinates on the axes
            gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), color='gray', linestyle='--', alpha=0.5)
            gl.top_labels = False    # Hides the labels on the top axis for a cleaner look
            gl.right_labels = False  # Hides the labels on the right axis
            
            plt.title(f"Arrival Times & Multiple Hits for {storm}")
            
            plt.show()
            
            
            ## Test plot showing SST over the storm period ##
            
            # 1. Define the path to the resampled SST file based on your other script
            sst_file_path = path.join(f"maxss\\storm-atlas\\ibtracs\\{region}\\{year}\\{storm}\\Resampled_for_fluxengine_MAXSS_ESACCI_SST.nc")
            
            # 2. Check if the SST file actually exists before trying to plot
            if os.path.exists(sst_file_path):
                # Load the SST NetCDF dataset
                sst_nc = nc.Dataset(sst_file_path)
                
                # 3. Extract the SST time series for the SAME grid cell we used for the wind plot
                # We reuse test_lat and test_lon here
                sst_ts = sst_nc.variables['sst'][:, test_lat, test_lon]
                
                # 4. Create the plot!
                plt.figure(figsize=(14, 6))
                
                # Plot the SST against the FULL timeline
                # We reuse plot_times which we already converted to standard datetimes
                plt.plot(plot_times, sst_ts, label='Sea Surface Temp (SST)', color='darkred', linewidth=1.5)
                
                # Shade the background using the arrival and departure datetimes directly
                if t_arrive_dt and t_depart_dt:
                    plt.axvspan(t_arrive_dt, t_depart_dt, color='grey', alpha=0.3, label='Active Storm Window')
                
                # Plot the timeline boundaries (reusing the exact same dt variables)
                if t_pre_start_dt: plt.axvline(t_pre_start_dt, color='blue', linestyle='--', linewidth=2, label='Pre-Storm Start (-15d)')
                if t_pre_end_dt: plt.axvline(t_pre_end_dt, color='cyan', linestyle='--', linewidth=2, label='Pre-Storm End (-2d)')
                if t_arrive_dt: plt.axvline(t_arrive_dt, color='red', linestyle='-', linewidth=2.5, label='First Arrival')
                if t_depart_dt: plt.axvline(t_depart_dt, color='purple', linestyle='-', linewidth=2.5, label='Final Departure')
                if t_post_end_dt: plt.axvline(t_post_end_dt, color='green', linestyle='--', linewidth=2, label='Post-Storm End (+40d)')
                
                # Formatting
                plt.title(f"Diagnostic SST Timeline for {storm} (Lat: {matched_lat}, Lon: {matched_lon})")
                plt.ylabel("Sea Surface Temperature (Kelvin)")
                plt.xlabel("Date")
                
                # Move legend outside the plot so it doesn't block the data
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                plt.grid(True, alpha=0.3)
                plt.tight_layout()
                plt.show()
                
                # Always close the NetCDF file when done!
                sst_nc.close()
            else:
                print(f"Skipping SST plot: Could not find SST file at {sst_file_path}")
            
            ## --- COMBINED DIAGNOSTIC TIMELINE: WIND & SST --- ##
            
            # 1. Define the figure and the left axis (for Wind Speed)
            fig, ax1 = plt.subplots(figsize=(14, 6))
            
            color1 = 'black'
            ax1.set_xlabel("Date")
            ax1.set_ylabel("Wind Speed (m/s)", color=color1)
            ax1.tick_params(axis='y', labelcolor=color1)
            
            # Plot the actual wind speed against the FULL timeline
            line1, = ax1.plot(plot_times, wind_ts, label='Wind Speed', color=color1, linewidth=1.5)
            
            # Shade the background using the arrival and departure datetimes directly
            if t_arrive_dt and t_depart_dt:
                ax1.axvspan(t_arrive_dt, t_depart_dt, color='grey', alpha=0.3, label='Active Storm Window')
    
            # Plot the timeline boundaries on ax1
            if t_pre_start_dt: ax1.axvline(t_pre_start_dt, color='blue', linestyle='--', linewidth=2, label='Pre-Storm Start (-15d)')
            if t_pre_end_dt: ax1.axvline(t_pre_end_dt, color='cyan', linestyle='--', linewidth=2, label='Pre-Storm End (-2d)')
            if t_arrive_dt: ax1.axvline(t_arrive_dt, color='red', linestyle='-', linewidth=2.5, label='First Arrival')
            if t_depart_dt: ax1.axvline(t_depart_dt, color='purple', linestyle='-', linewidth=2.5, label='Final Departure')
            if t_post_end_dt: ax1.axvline(t_post_end_dt, color='green', linestyle='--', linewidth=2, label='Post-Storm End (+40d)')
            
            # 2. Create the right axis (for SST) sharing the same X-axis
            ax2 = ax1.twinx()
            
            color2 = 'darkred'
            ax2.set_ylabel("Sea Surface Temperature (Kelvin)", color=color2)
            ax2.tick_params(axis='y', labelcolor=color2)
            
            # 3. Check for the SST file and extract the data for the exact same cell
            sst_file_path = path.join(f"maxss\\storm-atlas\\ibtracs\\{region}\\{year}\\{storm}\\Resampled_for_fluxengine_MAXSS_ESACCI_SST.nc")
            line2 = None
            if os.path.exists(sst_file_path):
                sst_nc = nc.Dataset(sst_file_path)
                sst_ts = sst_nc.variables['sst'][:, test_lat, test_lon]
                line2, = ax2.plot(plot_times, sst_ts, label='Sea Surface Temp (SST)', color=color2, linewidth=1.5)
                sst_nc.close()
            else:
                print(f"Warning: SST file not found for {storm}. Plotting wind only.")

            # 4. Combine legends from both axes so they don't overlap
            handles, labels = ax1.get_legend_handles_labels()
            if line2:
                handles.append(line2)
                labels.append(line2.get_label())
                
            # Move combined legend outside the plot
            ax1.legend(handles, labels, loc='center left', bbox_to_anchor=(1.1, 0.5))
            
            # Formatting
            plt.title(f"Diagnostic Timeline (Wind & SST) for {storm} (Lat: {matched_lat}, Lon: {matched_lon})")
            ax1.grid(True, alpha=0.3)
            
            # tight_layout ensures nothing gets cut off by the new legend position
            fig.tight_layout()
            plt.show()
            
            ## --- END COMBINED DIAGNOSTIC TIMELINE --- ##        
            
            
            # Remember to close your NetCDF file at the end of the loop to free up memory!
            winds_nc.close()
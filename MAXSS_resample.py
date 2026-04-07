# -*- coding: utf-8 -*-
"""

Script for looping through the MAXSS storm-Atlas and resampling data for fluxengine runs

@Original author: Richard Sims

This script has been updated in 2026 by Daniel Wilson

"""

#### all functions and scripts imported here

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
import datetime


def process_slice(valData, errData, outputRes=0.25):
    
    #Handle masked values explicitly to avoid masked element to nan error
    valData = valData.filled(fill_value=np.nan)
    errData = errData.filled(fill_value=np.nan)
    
    newGrid = np.full((wind_lat_dimension, wind_lon_dimension), np.nan, dtype=float);
    newGridCount = np.zeros((wind_lat_dimension, wind_lon_dimension), dtype=float);
    newGridErr = np.full((wind_lat_dimension, wind_lon_dimension), np.nan, dtype=float);
    
    #Supress warning for landmass/empty cells
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
    
        for ilat, lat in enumerate(np.arange(min_lat, max_lat, outputRes)):
            for ilon, lon in enumerate(np.arange(min_lon, max_lon, outputRes)):
                if iCoordMeshes[ilat, ilon] is not None:
                    newGrid[ilat, ilon] = np.nanmean(valData[iCoordMeshes[ilat, ilon]]);
                    #newGridCount[ilat, ilon] = np.nansum(countData[iCoordMeshes[ilat, ilon]]);
                    newGridCount[ilat, ilon] = len(iCoordMeshes[ilat, ilon][0]) * len(iCoordMeshes[ilat, ilon][0][0]);
                    newGridErr[ilat, ilon] = np.sqrt(np.nansum(pow(errData[iCoordMeshes[ilat, ilon]],2)));
    
    newGridErr[newGridCount!=0] = newGridErr[newGridCount!=0] / newGridCount[newGridCount!=0];
    
    return newGrid, newGridCount, newGridErr;


if __name__ == "__main__":

#### Get the path of the root directory where the data are stored. 
    #This will be user specific and can be changed depening on where data is stored
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
                if any(name in storm for name in [ "BONNIE", "COLIN", "MARIA", "ALEX" ]): #"RINA
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
                    
                    
                #### load L4 Wind data
                winds_nc = nc.Dataset(path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\MAXSS_HIST_TC_{3}_{1}_{4}_MAXSS_HIST_TC_L4.nc".format(region,year,storm,region_id,storm_id)));
                
                # 1. Extract the required variables from the NetCDF file
                spatial_mask = winds_nc.variables['data_spatial_mask'][:] 
                track_times = winds_nc.variables['__track_time'][:]
                track_time_units = winds_nc.variables['__track_time'].units # Needed for datetime conversion
                
                # 2. Grab the EXACT fill value from the metadata
                # This ensures your output matches the input perfectly
                fill_value = winds_nc.variables['__eo_eastward_wind']._FillValue
                
                # 3. Find out if the storm EVER passed over each grid cell
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
                
                ## Create 3D Analysis Period Mask ##
                
                # 1. Grab the full data time array 
                full_data_times = winds_nc.variables['time'][:]
                full_data_units = winds_nc.variables['time'].units
                
                # 2. Reshape for 3D broadcasting
                full_data_times_3d = full_data_times[:, np.newaxis, np.newaxis]
                
                # 3. Re-create the 3D mask using the FULL timeline
                # This compares every hour of the data to your spatial 2D arrival/departure maps
                analysis_period_mask_3d = (full_data_times_3d >= pre_storm_ref_period_start_numeric) & \
                                          (full_data_times_3d <= post_storm_analysis_end_numeric)
                
                # Ensure this mask only keeps areas hit by the storm
                analysis_period_mask_3d = np.logical_and(analysis_period_mask_3d, ever_in_storm_mask)
                
                #print(f"New Mask Shape: {analysis_period_mask_3d.shape}") 
                
                ## Calculating the pre-storm reference ##
                
                # 1. This compares the full timeline against the 15-day and 2-day offsets
                pre_storm_mask_3d = (full_data_times_3d >= pre_storm_ref_period_start_numeric) & \
                                    (full_data_times_3d <= pre_storm_ref_period_end_numeric)
                                    
                # 2. Only keep values for pixels the storm actually hit
                pre_storm_mask_3d = np.logical_and(pre_storm_mask_3d, ever_in_storm_mask)
                
                ## Save the storm timings and analysis period mask to a netcdf file ##
                
                # 1. Load in required variables from wind netcdf
                raw_eastward = winds_nc.variables['__eo_eastward_wind'][:]
                raw_northward = winds_nc.variables['__eo_northward_wind'][:]
                
                # 2. Convert the mask to actual NaNs and THEN change to float32
                # .filled(np.nan) ensures that the hidden fill values become 'Not a Number'
                # so they are ignored by your later np.nanmedian() calls.
                wind_eastward = raw_eastward.filled(np.nan).astype('float32')
                wind_northward = raw_northward.filled(np.nan).astype('float32')
                
                # 2. get wind data dimensions
                wind_time_dimension=len(wind_eastward)
                wind_lat_dimension=len(wind_eastward[0])
                wind_lon_dimension=len(wind_eastward[0][0])
                
                # 3. define lat,long,time variables - these are used to build new matrixes
                wind_lat = winds_nc.variables['lat'][:]
                wind_lon = winds_nc.variables['lon'][:]
                wind_time = winds_nc.variables['time'][:]
                wind_dates = num2date(wind_time, winds_nc.variables['time'].units)
                 
                # 4. min and max of wind grid lats and lons   
                min_lat=wind_lat[0]
                max_lat=wind_lat[0-1]                    
                min_lon=wind_lon[0]
                max_lon=wind_lon[0-1] 
                
                #5. set up netcdf
                processedFilePath = (path.join(f"maxss\\storm-atlas\\ibtracs\\{region}\\{year}\\{storm}\\Resampled_for_fluxengine_storm_timings_with_masks.nc"))
                ncout = Dataset(processedFilePath, 'w')
                
                # 6. Define Dimensions
                # We now include 'time' which matches your 475-step data timeline
                ncout.createDimension("time", len(full_data_times))
                ncout.createDimension("lat", wind_lat_dimension)
                ncout.createDimension("lon", wind_lon_dimension)
                
                # 7. Create Dimension Variables
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

                # 8. Create the 3D Mask Variable
                # Use 'int8' (or 'byte') to save space since it's just 0s and 1s
                var_mask = ncout.createVariable("analysis_mask", "i1", ("time", "lat", "lon"), 
                                                zlib=True, complevel=4, shuffle=True)
                var_mask.long_name = "3D Analysis Period Mask (-15 to +40 days)"
                var_mask.description = "1 = Valid window for that pixel, 0 = Outside window"
                
                # 10. We convert the boolean mask to integers (0 and 1) for saving
                var_mask[:] = analysis_period_mask_3d.astype(np.int8)
                
                # 11. Create the 3D Mask Variable
                # Use 'int8' (or 'byte') to save space since it's just 0s and 1s
                pre_storm_ref_mask_var = ncout.createVariable("pre_storm_ref_mask", "i1", ("time", "lat", "lon"), 
                                                zlib=True, complevel=4, shuffle=True)
                pre_storm_ref_mask_var.long_name = "3D Pre Storm reference Period Mask (-15 to -2 days)"
                pre_storm_ref_mask_var.description = "1 = Valid window for that pixel, 0 = Outside window"
                
                # 12. We convert the boolean mask to integers (0 and 1) for saving
                pre_storm_ref_mask_var[:] = pre_storm_mask_3d.astype(np.int8)

                # 13. Create the 2D Timing Variables (Same as your original code)
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
                
                ## MAXSS L4 Wind data (wind speed and moment 2) ##
                
                # calculate wind speed and wind moment2 from east and west components
                wind_speed = np.hypot(wind_eastward, wind_northward)
                wind_moment2 = wind_speed**2
                
                # Any pixel where the mask is 0 (outside the -15 to +40 window) 
                # gets converted to Not-a-Number (NaN)
                wind_speed[analysis_period_mask_3d == 0] = np.nan
                wind_moment2[analysis_period_mask_3d == 0] = np.nan
                
                # Check if data is masked array or standard narray and standardise by replacing nans with MAXSS fill value.
                if hasattr(wind_speed, 'filled'):
                    wind_speed = wind_speed.filled(fill_value=fill_value).astype('float32')
                else:
                    wind_speed = np.nan_to_num(wind_speed, nan=fill_value).astype('float32')
                
                if hasattr(wind_moment2, 'filled'):
                    wind_moment2 = wind_moment2.filled(fill_value=fill_value).astype('float32')
                else:
                    wind_moment2 = np.nan_to_num(wind_moment2, nan=fill_value).astype('float32')
                
                # save wind speed to netCDF
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_MAXSS_L4_windspeed.nc".format(region,year,storm)));

                ncout = Dataset(processedFilePath, 'w');
                
                # provide dimensions
                ncout.createDimension("lat", wind_lat_dimension);
                ncout.createDimension("lon", wind_lon_dimension);
                ncout.createDimension("time", wind_time_dimension);
                
                # dimension variables
                var = ncout.createVariable("lat", float, ("lat",));
                var.units = "lat (degrees North)";
                var[:] = wind_lat;
                
                var = ncout.createVariable("lon", float, ("lon",));
                var.units = "lon (degrees East)";
                var[:] = wind_lon;
                
                var = ncout.createVariable("time", int, ("time",));
                var.long_name = "Time";
                var.units = "seconds since 1981-01-01";
                var[:] = wind_time
                
                #data variables
                var = ncout.createVariable("windspeed", "f4", ("time","lat", "lon"), 
                                           zlib=True, complevel=1, shuffle=True, fill_value=fill_value,
                                           chunksizes=(1, wind_lat_dimension, wind_lon_dimension));
                var.units = "m s-1";
                var.long_name = "Hourly wind speed from MAXSS on a 0.25X0.25 degree gridspatial resolution";
                var[:] = wind_speed;
                
                #data variables
                var = ncout.createVariable("second_moment_wind", "f4", ("time","lat", "lon"), 
                                           zlib=True, complevel=1, shuffle=True, fill_value=fill_value,
                                           chunksizes=(1, wind_lat_dimension, wind_lon_dimension));
                var.units = "m2 s-2";
                var.long_name = "Second moment of wind speed from MAXSS on a 0.25X0.25 degree gridspatial resolution";
                var[:] = wind_moment2;
                
                ncout.close();  
                
                ## MAXSS L4 Wind: Median pre storm conditions over the 13 day pre-storm period ##
                               
                # 1. Create a copy of the data as floats so it can safely hold 'NaN' values
                pre_storm_data = wind_speed.astype('float32')
                pre_storm_moment_data = wind_moment2.astype('float32') 
                
                # 2. Apply the pre-storm mask! 
                # "Anywhere the pre_storm_mask_3d is False, change the data to NaN"
                pre_storm_data[pre_storm_mask_3d == 0] = np.nan
                pre_storm_moment_data[pre_storm_mask_3d == 0] = np.nan
                
                # 3. Calculate the median along the time axis (axis=0)
                # np.nanmedian ignores the NaNs, so it ONLY calculates the median of the 
                # valid data points inside your -15 to -2 day window!
                # (We use warnings.catch_warnings to suppress the warning NumPy throws when 
                # it encounters a grid cell that is 100% NaNs, like landmasses or areas the storm missed)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    pre_storm_median_2d = np.nanmedian(pre_storm_data, axis=0)
                    pre_storm_moment_median_2d = np.nanmedian(pre_storm_moment_data, axis=0)
                
                print(f"Calculated 2D pre-storm median data for {storm}.")
                
                # 4. Broadcast the 2D arrays into 3D hourly format for Fluxengine!
                # np.tile instantly copies the 2D map across the entire time dimension
                wind_speed_prestormref_3d = np.tile(pre_storm_median_2d, (wind_time_dimension, 1, 1))
                wind_moment2_prestormref_3d = np.tile(pre_storm_moment_median_2d, (wind_time_dimension, 1, 1))
                
                # 5. Apply the 3D analysis period mask to remove data that will not be used in Flux calculations
                wind_speed_prestormref_3d[analysis_period_mask_3d == 0] = np.nan
                wind_moment2_prestormref_3d[analysis_period_mask_3d == 0] = np.nan
                
                # 6. Replace NaNs with source fill value (converts both land and areas we just masked)
                # This ensures that when you TILE them, the 3D array is already "clean"
                wind_speed_prestormref_3d = np.nan_to_num(wind_speed_prestormref_3d, nan=fill_value).astype('float32')
                wind_moment2_prestormref_3d = np.nan_to_num(wind_moment2_prestormref_3d, nan=fill_value).astype('float32')
                
                # 7. Save wind pre-storm output into a netCDF 
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_MAXSS_L4_windspeed_pre_storm_reference.nc".format(region,year,storm)))

                ncout = Dataset(processedFilePath, 'w')
                    
                # Create dataset and provide dimensions
                ncout.createDimension("lat", wind_lat_dimension)
                ncout.createDimension("lon", wind_lon_dimension)
                ncout.createDimension("time", wind_time_dimension)
                
                # Dimension variables
                var = ncout.createVariable("lat", float, ("lat",))
                var.units = "lat (degrees North)"
                var[:] = wind_lat
                
                var = ncout.createVariable("lon", float, ("lon",))
                var.units = "lon (degrees East)"
                var[:] = wind_lon
                
                var = ncout.createVariable("time", int, ("time",))
                var.long_name = "Time"
                var.units = "seconds since 1981-01-01"
                var[:] = wind_time
                
                # Data variables
                var = ncout.createVariable("windspeed", "f4", ("time","lat", "lon"), 
                                           zlib=True, complevel=1, shuffle=True, fill_value=fill_value) #chunksizes=(1, wind_lat_dimension, wind_lon_dimension)
                                           
                var.units = "m s-1"
                var.long_name = "Dynamic Pre-storm (-15 to -2 days) median of hourly wind speed"
                var[:] = wind_speed_prestormref_3d
                
                var = ncout.createVariable("second_moment_wind", "f4", ("time","lat", "lon"), 
                                           zlib=True, complevel=1, shuffle=True, fill_value=fill_value)#;chunksizes=(1, wind_lat_dimension, wind_lon_dimension)
                var.units = "m2 s-2";
                var.long_name = "Dynamic Pre-storm (-15 to -2 days) median of Second moment of wind speed from MAXSS on a 0.25X0.25 degree gridspatial resolution";
                var[:] = wind_moment2_prestormref_3d;
                
                ncout.close()
                print(f"Successfully saved hourly Fluxengine reference map for {storm}.")
                
                del wind_eastward, wind_northward  # Large 3D input arrays
                del wind_speed, wind_moment2       # Large 3D calculated arrays
                del pre_storm_data, pre_storm_moment_data # Large 3D masked copies
                del wind_speed_prestormref_3d, wind_moment2_prestormref_3d # Tiled 3D outputs
                
                gc.collect()
                
                print("Wind regridded for Storm = "+storm)
                
                
                ## MAXSS ESACCI SST data ##
                # load data
                sst_nc = nc.Dataset(path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\MAXSS_HIST_TC_{3}_{1}_{4}_OSTIA-C3S-L4-GLOB-v2.1.nc".format(region,year,storm,region_id,storm_id)));
                #sst = sst_nc.variables['__eo_analysed_sst'][:]
                sst_lat = sst_nc.variables['lat'][:]
                sst_lon = sst_nc.variables['lon'][:]
                sst_time = sst_nc.variables['time'][:]
                
                # EXTRACT THE SST SPECIFIC FILL VALUE
                sst_fill_value = sst_nc.variables['__eo_analysed_sst']._FillValue
                
                # resample SST to wind grid
                iCoordMeshes = None; #Initially this is None but will be calculated exactly once. It's a long calculation so doesn't want to be repeated.
                outputRes = 0.25;
                # Calculate binning information
                #Only do this once because it's computationally expensive but the same for all time steps
                if iCoordMeshes is None: 
                    CCILats = sst_nc.variables['lat'][:]
                    CCILons = sst_nc.variables['lon'][:]
                    #print("Calculating grid cell mapping...");
                    iCoordMeshes = np.full((wind_lat_dimension,wind_lon_dimension), None, dtype=object);
                    
                    for ilat, lat in enumerate(np.arange(min_lat,max_lat , outputRes)):
                        #print("Grid cell mapping for latitude", lat);
                        for ilon, lon in enumerate(np.arange(min_lon,max_lon, outputRes)):
                            wlat = np.where((CCILats >= lat) & (CCILats < (lat+outputRes)));
                            wlon = np.where((CCILons >= lon) & (CCILons< (lon+outputRes)));
                            
                            if (len(wlat[0]) > 0) & (len(wlon[0]) > 0):
                                iCoordMeshes[ilat, ilon] = np.meshgrid(wlat[0], wlon[0]);
                                    
                # loop through timesteps 
                timesteps_sst=len(sst_time)
                del ilat, ilon, CCILats,CCILons, lat , lon
                
                # Store data for each day of the month
                sst_regrid_Vals = np.empty((timesteps_sst, wind_lat_dimension,wind_lon_dimension), dtype=float);
                sst_regrid_ValsErr = np.empty((timesteps_sst, wind_lat_dimension,wind_lon_dimension), dtype=float);
                sst_regrid_ValsCounts = np.empty((timesteps_sst, wind_lat_dimension,wind_lon_dimension), dtype=float);
                               
                for sst_timesteps in range(0, timesteps_sst): 
                    #print(sst_timesteps)                               
                    sst_time_slice=sst_nc.variables['__eo_analysed_sst'][sst_timesteps,:,:]                         
                    sst_uncertainty_slice=sst_nc.variables['__eo_analysed_sst_uncertainty'][sst_timesteps,:,:]                         
                    newVals, newCountCount, newValsErr = process_slice(sst_time_slice, sst_uncertainty_slice);
                
                    sst_regrid_Vals[sst_timesteps,:,:] = newVals;
                    sst_regrid_ValsErr[sst_timesteps,:,:] = newValsErr;
                    sst_regrid_ValsCounts[sst_timesteps,:,:] = newCountCount;
                
                #plt.pcolor(sst_regrid_Vals[1,:,:])
                             
                sst_on_wind_grid = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                
                # get the data of each timestep in wind data
                sst_time = sst_nc.variables['time'][:]
                sst_dates = num2date(sst_time, sst_nc.variables['time'].units)
                
                # loop through the wind timestamps, extract the month and use that to pick which SST data to use.
                for wind_step in range(0, wind_time_dimension):
                
                    #now find the index of the CLOSEST value and use that
                    abs_deltas_from_target_date = np.absolute(sst_dates - wind_dates[wind_step])
                    index_of_min_delta_from_target_date = np.argmin(abs_deltas_from_target_date)
                    #closest_date = sst_dates[index_of_min_delta_from_target_date]

                    sst_on_wind_grid[wind_step,:,:]=sst_regrid_Vals[index_of_min_delta_from_target_date,:,:]
                    
                # Apply the 3D analysis period mask (-15 to +40 days)
                #(This mask already inherently includes the 'ever_in_storm' boundaries)
                sst_on_wind_grid[analysis_period_mask_3d == 0] = np.nan
                
                #convert data to fill values rather than nan
                sst_on_wind_grid = np.nan_to_num(sst_on_wind_grid, nan=sst_fill_value).astype('float32')
                
                # save SST output into a netCDF 
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_MAXSS_ESACCI_SST.nc".format(region,year,storm)));

                ncout = Dataset(processedFilePath, 'w');
                    
                # create dataset and provide dimensions
                ncout.createDimension("lat", wind_lat_dimension);
                ncout.createDimension("lon", wind_lon_dimension);
                ncout.createDimension("time", wind_time_dimension);
                
                # dimension variables
                var = ncout.createVariable("lat", "f4", ("lat",));
                var.units = "lat (degrees North)";
                var[:] = wind_lat;
                
                var = ncout.createVariable("lon", "f4", ("lon",));
                var.units = "lon (degrees East)";
                var[:] = wind_lon;
                
                var = ncout.createVariable("time", int, ("time",));
                var.long_name = "Time";
                var.units = "seconds since 1981-01-01";
                var[:] = wind_time
                
                # data variables
                var = ncout.createVariable("sst", "f4", ("time","lat", "lon"), 
                                           zlib=True, complevel=1, shuffle=True, fill_value=sst_fill_value);
                var.units = "Degrees Kelvin";
                var.long_name = "Daily ESACCI sea surface temperature resampled to a 0.25X0.25 degree spatial and hourly temporal resolutionn";
                var[:] = sst_on_wind_grid;
                
                ncout.close();   
                
                
                ## MAXSS ESACCI SST Median pre storm conditions over the 13 day pre-storm period ##
                           
                # 1. Create a float copy for safe NaN handling
                pre_storm_sst = sst_on_wind_grid.copy()
                
                # 2. Apply the dynamic 3D mask you already calculated!
                pre_storm_sst[pre_storm_mask_3d == 0] = np.nan
                
                # 3. Calculate the median/mean ONLY inside the pre-storm window
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    # Using nanmedian to match wind methodology
                    SST_prestormref_2d = np.nanmedian(pre_storm_sst, axis=0)
                    
                # 4. Tile the data into 3D
                sst_on_wind_grid_prestormref = np.tile(SST_prestormref_2d, (wind_time_dimension, 1, 1))
                
                # 5. Apply the 3D analysis period mask
                sst_on_wind_grid_prestormref[analysis_period_mask_3d == 0] = np.nan
                
                # 6. Clean up NaNs with the correct fill value
                sst_on_wind_grid_prestormref = np.nan_to_num(sst_on_wind_grid_prestormref, nan=sst_fill_value).astype('float32')
                
                # 7. Save SST pre storm median into a netCDF 
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_MAXSS_ESACCI_SST_pre_storm_reference.nc".format(region,year,storm)));

                ncout = Dataset(processedFilePath, 'w');
                    
                # create dataset and provide dimensions
                ncout.createDimension("lat", wind_lat_dimension);
                ncout.createDimension("lon", wind_lon_dimension);
                ncout.createDimension("time", wind_time_dimension);
                
                # dimension variables
                var = ncout.createVariable("lat", "f4", ("lat",));
                var.units = "lat (degrees North)";
                var[:] = wind_lat;
                
                var = ncout.createVariable("lon", "f4", ("lon",));
                var.units = "lon (degrees East)";
                var[:] = wind_lon;
                
                var = ncout.createVariable("time", int, ("time",));
                var.long_name = "Time";
                var.units = "seconds since 1981-01-01";
                var[:] = wind_time
                
                # data variables
                var = ncout.createVariable("sst", "f4", ("time","lat", "lon"), 
                                           zlib=True, complevel=1, shuffle=True, fill_value=sst_fill_value);
                var.units = "Degrees Kelvin";
                var.long_name = "Dynamic Pre-storm (-15 to -2 days) median of daily ESACCI sea surface temperature resampled to a 0.25X0.25 degree spatial and hourly temporal resolution";
                var[:] = sst_on_wind_grid_prestormref;
                
                ncout.close();   
                
                # 8. Safely clean up un-used variables
                # Combine all variables into a single list and check if they exist before deleting
                vars_to_delete = [
                    'sst_on_wind_grid', 'sst_regrid_Vals', 'pre_storm_sst', 
                    'sst_on_wind_grid_prestormref', 'SST_prestormref_2d', 
                    'sst_dates', 'sst_lat', 'sst_lon', 'sst_nc', 
                    'newVals', 'sst_regrid_ValsErr', 'sst_regrid_ValsCounts', 
                    'newCountCount', 'newValsErr', 'sst_time', 'sst_time_slice', 
                    'sst_uncertainty_slice', 'timesteps_sst', 'abs_deltas_from_target_date', 
                    'index_of_min_delta_from_target_date']
                
                for var_name in vars_to_delete:
                    if var_name in locals():
                        del locals()[var_name]
                
                gc.collect()
                print("SST regridded for Storm = "+storm)
         

                ## MAXSS ESACCI SSS data ##
                # 1. load data                
                sss_nc = nc.Dataset(path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\MAXSS_HIST_TC_{3}_{1}_{4}_ESACCI-SEASURFACESALINITY-L4-SSS-MERGED_OI_7DAY_RUNNINGMEAN_DAILY_25km.nc".format(region,year,storm,region_id,storm_id)));
                #sss = sss_nc.variables['__eo_sss'][:]
                sss_lat = sss_nc.variables['lat'][:]
                sss_lon = sss_nc.variables['lon'][:]
                sss_time = sss_nc.variables['time'][:]
                
                    #### resample SSS to wind grid
                iCoordMeshes = None; #Initially this is None but will be calculated exactly once. It's a long calculation so doesn't want to be repeated.
                outputRes = 0.25;
                #Only do this once because it's computationally expensive but the same for all time steps
                    #### Calculate binning information
                if iCoordMeshes is None: 
                    CCILats = sss_nc.variables['lat'][:]
                    CCILons = sss_nc.variables['lon'][:]
                    #print("Calculating grid cell mapping...");
                    iCoordMeshes = np.full((wind_time_dimension, wind_lon_dimension), None, dtype=object);
                
                    for ilat, lat in enumerate(np.arange(min_lat,max_lat , outputRes)):
                        #print("Grid cell mapping for latitude", lat);
                        for ilon, lon in enumerate(np.arange(min_lon,max_lon, outputRes)):
                            wlat = np.where((CCILats >= lat) & (CCILats < (lat+outputRes)));
                            wlon = np.where((CCILons >= lon) & (CCILons< (lon+outputRes)));
                            
                            if (len(wlat[0]) > 0) & (len(wlon[0]) > 0):
                                iCoordMeshes[ilat, ilon] = np.meshgrid(wlat[0], wlon[0]);
                                    
                    # loop through timesteps 
                timesteps_sss=len(sss_time)
                
                    #### Store data for each week of the month
                sss_regrid_Vals = np.empty((timesteps_sss, wind_lat_dimension,wind_lon_dimension), dtype=float);
                #sss_regrid_ValsErr = np.empty((timesteps_sss, wind_lat_dimension,wind_lon_dimension), dtype=float);
                #sss_regrid_ValsCounts = np.empty((timesteps_sss, wind_lat_dimension,wind_lon_dimension), dtype=float);
                               
                for sss_timesteps in range(0, timesteps_sss): 
                    #print(sss_timesteps)                               
                    sss_time_slice=sss_nc.variables['__eo_sss'][sss_timesteps,:,:]                         
                    sss_uncertainty_slice=sss_nc.variables['__eo_sss_random_error'][sss_timesteps,:,:]                         
                    newVals, newCountCount, newValsErr = process_slice(sss_time_slice, sss_uncertainty_slice);
                
                    #regridding results in some empty tracks
                    x_coord_range = wind_lon.tolist()
                    y_coord_range = wind_lat.tolist()
                    xy_coord = list(itertools.product(x_coord_range, y_coord_range))  
                    #this is first input into interpolation function
                    sample_df = pd.DataFrame()
                    sample_df['X'] = [xy[0] for xy in xy_coord]
                    sample_df['Y'] = [xy[1] for xy in xy_coord]
                    x_sss=np.array([sample_df["X"]])
                    x_sss=x_sss.transpose()
                    x_sss = np.squeeze(x_sss)
                    y=np.array([sample_df["Y"]])
                    y=y.transpose()
                    y = np.squeeze(y)
                    values = newVals.flatten(order='F')
                    #this is second input into interpolation function
                    sample_df1 = pd.DataFrame()
                    sample_df1['value'] = values
                    z=np.array([sample_df1['value']])
                    z=z.transpose()
                    z = np.squeeze(z)
                    index=np.argwhere(np.isnan(z))
                    x_sss=np.delete(x_sss, index)
                    y=np.delete(y, index)
                    z=np.delete(z, index)
                    grd_lon = wind_lon
                    grd_lat = wind_lat
                    grd_lons, grd_lats = np.meshgrid(grd_lon, grd_lat)
                    SSS_nearest = griddata((x_sss, y), z, (grd_lons, grd_lats), method='nearest')
                
                    
                    #now remove the stuff happening at the edges by making 
                    #this is from SST grid boolean True
                    where_sst_nan=np.isnan(sst_regrid_Vals[0,:,:])
                    xcv=np.where(where_sst_nan)
                    p=xcv[0]
                    pp=xcv[1]
                    SSS_nearest[p,pp]=float("nan")
                
                    #SSS_nearest=np.flipud(SSS_nearest)
                    sss_regrid_Vals[sss_timesteps,:,:] = SSS_nearest;
                    #sss_regrid_ValsErr[sss_timesteps,:,:] = newValsErr;
                    #sss_regrid_ValsCounts[sss_timesteps,:,:] = newCountCount;
                
                
                    #### get the data of each timestep in wind data
                sss_on_wind_grid = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                sss_time = sss_nc.variables['time'][:]
                sss_dates = num2date(sss_time, sss_nc.variables['time'].units)
                
                    ####loop through the wind timestamps
                    #extract the dates and find the closest SST data.
                for wind_step in range(0, wind_time_dimension):
                
                    #now find the index of the CLOSEST value and use that
                    abs_deltas_from_target_date = np.absolute(sss_dates - wind_dates[wind_step])
                    index_of_min_delta_from_target_date = np.argmin(abs_deltas_from_target_date)
                    #closest_date = sss_dates[index_of_min_delta_from_target_date]
                
                    
                    sss_on_wind_grid[wind_step,:,:]=sss_regrid_Vals[index_of_min_delta_from_target_date,:,:]
                    
                
                    #### save SSS output into a netCDF 
                
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_MAXSS_ESACCI_SSS.nc".format(region,year,storm)));
                ncout = Dataset(processedFilePath, 'w');
                    #### create dataset and provide dimensions
                
                ncout.createDimension("lat", wind_lat_dimension);
                ncout.createDimension("lon", wind_lon_dimension);
                ncout.createDimension("time", wind_time_dimension);
                
                #dimension variables
                var = ncout.createVariable("lat", float, ("lat",));
                var.units = "lat (degrees North)";
                var[:] = wind_lat;
                
                var = ncout.createVariable("lon", float, ("lon",));
                var.units = "lon (degrees East)";
                var[:] = wind_lon;
                
                var = ncout.createVariable("time", int, ("time",));
                var.long_name = "Time";
                var.units = "seconds since 1981-01-01";
                var[:] = wind_time
                
                #data variables
                var = ncout.createVariable("sss", float, ("time","lat", "lon"), 
                                           zlib=True, complevel=1, shuffle=True, chunksizes=(1, wind_lat_dimension, wind_lon_dimension));
                var.units = "PSU";
                var.long_name = "Weekly ESACCI sea surface salinity resampled to a 0.25X0.25 degree spatial and hourly temporal resolution";
                var[:] = sss_on_wind_grid;
                
                ncout.close();   
                
                #### MAXSS ESACCI SSS data pre_storm_reference

                #pre storm is first 15 days,so 15 days * hourly resolution
                SSS_prestormref=np.nanmean(sss_on_wind_grid[0:(15*24):1,:,:],axis =(0))
                #create empty grid
                sss_on_wind_grid_prestormref = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                #set all the values equal to the pre storm mean values
                for wind_step in range(0, wind_time_dimension):
                    sss_on_wind_grid_prestormref[wind_step,:,:]=SSS_prestormref[:,:]
                
                    #### save SSS output into a netCDF 
                
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_MAXSS_ESACCI_SSS_pre_storm_reference.nc".format(region,year,storm)));
                ncout = Dataset(processedFilePath, 'w');
                    #### create dataset and provide dimensions
                
                ncout.createDimension("lat", wind_lat_dimension);
                ncout.createDimension("lon", wind_lon_dimension);
                ncout.createDimension("time", wind_time_dimension);
                
                #dimension variables
                var = ncout.createVariable("lat", float, ("lat",));
                var.units = "lat (degrees North)";
                var[:] = wind_lat;
                
                var = ncout.createVariable("lon", float, ("lon",));
                var.units = "lon (degrees East)";
                var[:] = wind_lon;
                
                var = ncout.createVariable("time", int, ("time",));
                var.long_name = "Time";
                var.units = "seconds since 1981-01-01";
                var[:] = wind_time
                
                #data variables
                var = ncout.createVariable("sss", float, ("time","lat", "lon"), 
                                           zlib=True, complevel=1, shuffle=True, chunksizes=(1, wind_lat_dimension, wind_lon_dimension));
                var.units = "PSU";
                var.long_name = "Pre storm (15 days) mean of ESACCI weekly sea surface salinity resampled to a 0.25X0.25 degree spatial and hourly temporal resolution";
                var[:] = sss_on_wind_grid_prestormref;
                
                ncout.close();   
                
                    #### delete all the variables used during import and saving of sss      
                del SSS_prestormref,ncout, sss_dates, sss_lat,sss_lon,sss_nc,sss_on_wind_grid,sss_on_wind_grid_prestormref
                del newVals, sss_regrid_Vals,sss_time,sss_time_slice,sss_timesteps,sss_uncertainty_slice,SSS_nearest
                del processedFilePath,timesteps_sss,abs_deltas_from_target_date,iCoordMeshes,index_of_min_delta_from_target_date
                del where_sst_nan,values, newValsErr,newCountCount,p, pp, xcv, y, x_sss,z,xy_coord,grd_lat,grd_lon,grd_lats,grd_lons
                del sst_regrid_Vals,sample_df,sample_df1
                print("SSS regridded for Storm = "+storm)
          
                
                  
                #### MAXSS ERA5 precipitation 
                #Note this resample code assumes precipitation is on the same grid spatial and temporal grid as the wind
                #### load data
                precip_nc = nc.Dataset(path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\MAXSS_HIST_TC_{3}_{1}_{4}_ERA5_Precipitation.nc".format(region,year,storm,region_id,storm_id)))
                
                # Load data into memory
                p_lat = precip_nc.variables['lat'][:]
                p_lon = precip_nc.variables['lon'][:]
                p_time = precip_nc.variables['time'][:]
                
                # Load data and immediately replace mask with NaN to avoid issues with conversion
                p_data_raw = precip_nc.variables['__eo_tp'][:] # Shape: (time, lat, lon)
                p_data = np.ma.filled(p_data_raw, fill_value=np.nan)
                
                #### 1. Handle Latitude Orientation
                # If precip lat is [20, 19.75, 19.5] (descending) and wind lat is [19.5, 19.75, 20] (ascending)
                if p_lat[0] > p_lat[-1]:
                    p_lat = p_lat[::-1]
                    p_data = p_data[:, ::-1, :]
                    
                # 2. Check Shape/Dimensions 
                # Compare the lengths of the latitude and longitude arrays
                if (wind_lat.shape != p_lat.shape) or (wind_lon.shape != p_lon.shape):
                    print(f"\n!!! FATAL ERROR: Grid Dimension Mismatch for storm {storm} !!!")
                    print(f"  > Wind Grid Shape:   {wind_lat.shape} lat x {wind_lon.shape} lon")
                    print(f"  > Precip Grid Shape: {p_lat.shape} lat x {p_lon.shape} lon")
                    print("  > ACTION: Stopping script. You must regrid the inputs to match.")
                    sys.exit(1) # stops the script immediately
            
                # Check Coordinate Values
                # Use np.allclose to allow for tiny floating-point differences (e.g. 19.0000001 vs 19.0)
                # atol=1e-4 roughly equals 11 meters, which is strict enough for 0.25 deg grids
                lat_match = np.allclose(wind_lat, p_lat, atol=1e-4)
                lon_match = np.allclose(wind_lon, p_lon, atol=1e-4)
            
                if not lat_match or not lon_match:
                    print(f"\n!!! FATAL ERROR: Grid Coordinate Mismatch for storm {storm} !!!")
                    
                    if not lat_match:
                        print("  > Latitude arrays do not match.")
                        print(f"    Wind Lat Start/End:   {wind_lat[0]:.4f} / {wind_lat[-1]:.4f}")
                        print(f"    Precip Lat Start/End: {p_lat[0]:.4f} / {p_lat[-1]:.4f}")
                        
                        # Helpful Hint: Check if one is just flipped (Ascending vs Descending)
                        if np.allclose(wind_lat, p_lat[::-1], atol=1e-4):
                            print("    -> DIAGNOSIS: Latitudes are inverted (e.g., 90->-90 vs -90->90).")
                            print("    -> FIX: Enable the latitude flip code block in your script.")
            
                    if not lon_match:
                        print("  > Longitude arrays do not match.")
                        print(f"    Wind Lon Start/End:   {wind_lon[0]:.4f} / {wind_lon[-1]:.4f}")
                        print(f"    Precip Lon Start/End: {p_lon[0]:.4f} / {p_lon[-1]:.4f}")
                        
                        # Helpful Hint: Check for 0-360 vs -180-180 convention
                        if np.any(p_lon > 180) and np.all(wind_lon <= 180):
                            print("    -> DIAGNOSIS: Longitude conventions differ (0-360 vs -180-180).")
                    
                    print("  > ACTION: Stopping script to prevent physical errors.")
                    sys.exit(1)
            
                print(f"SUCCESS: Grids for {storm} are aligned.")
                
                #### 3. Vectorized Unit Conversion
                # ERA5 is 'm' per hour, FluxEngine needs 'mm/day'. 1m = 1000mm. Since it's hourly, * 24 for daily rate.
                unit_conversion = 1000 * 24
                p_data = p_data * unit_conversion

                #### 4. Temporal Alignment (Vectorized)
                precip_dates = num2date(p_time, precip_nc.variables['time'].units)
                precip_on_wind_grid = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=np.float32)
                
                for wind_step in range(wind_time_dimension):
                    # Find closest time index
                    abs_deltas = np.absolute(precip_dates - wind_dates[wind_step])
                    idx = np.argmin(abs_deltas)
                    
                    # Check if spatial dimensions match exactly. If so, just copy.
                    # If not, you can use simple slicing or scipy.interpolate.interp2d
                    precip_on_wind_grid[wind_step, :, :] = p_data[idx, :, :]
                
                print(f"Precipitation processed via vectorization for {storm}")
                          
                ####  5. Save precipitation pre storm output into a netCDF 
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_MAXSS_ERA5_precipitation.nc".format(region,year,storm)));
                ncout = Dataset(processedFilePath, 'w');
                    
                #### create dataset and provide dimensions
                ncout.createDimension("lat", wind_lat_dimension);
                ncout.createDimension("lon", wind_lon_dimension);
                ncout.createDimension("time", wind_time_dimension);
                
                #dimension variables
                var = ncout.createVariable("lat", float, ("lat",));
                var.units = "lat (degrees North)";
                var[:] = wind_lat;
                
                var = ncout.createVariable("lon", float, ("lon",));
                var.units = "lon (degrees East)";
                var[:] = wind_lon;
                
                var = ncout.createVariable("time", int, ("time",));
                var.long_name = "Time";
                var.units = "seconds since 1981-01-01";
                var[:] = wind_time
                
                #data variables
                var = ncout.createVariable("precipitation", float, ("time","lat", "lon"), 
                                           zlib=True, complevel=1, shuffle=True, chunksizes=(1, wind_lat_dimension, wind_lon_dimension));
                var.units = "mm day-1";
                var.long_name = "ERA5 hourly precipitation on a 0.25X0.25 degree spatial with hourly temporal resolution";
                var[:] = precip_on_wind_grid;
                
                ncout.close();   
                
                #### MAXSS ESACCI precipitation data pre_storm_reference
                #pre storm is first 15 days,so 15 days * hourly resolution
                precip_prestormref=np.nanmean(precip_on_wind_grid[0:(15*24):1,:,:],axis =(0))
                #create empty grid
                precip_on_wind_grid_prestormref = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                #set all the values equal to the pre storm mean values
                for wind_step in range(0, wind_time_dimension):
                    precip_on_wind_grid_prestormref[wind_step,:,:]=precip_prestormref[:,:]
                
                
                #### save precip output into a netCDF 
                
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_MAXSS_ERA5_precipitation_pre_storm_reference.nc".format(region,year,storm)));
                ncout = Dataset(processedFilePath, 'w');
                #### create dataset and provide dimensions
                
                ncout.createDimension("lat", wind_lat_dimension);
                ncout.createDimension("lon", wind_lon_dimension);
                ncout.createDimension("time", wind_time_dimension);
                
                #dimension variables
                var = ncout.createVariable("lat", float, ("lat",));
                var.units = "lat (degrees North)";
                var[:] = wind_lat;
                
                var = ncout.createVariable("lon", float, ("lon",));
                var.units = "lon (degrees East)";
                var[:] = wind_lon;
                
                var = ncout.createVariable("time", int, ("time",));
                var.long_name = "Time";
                var.units = "seconds since 1981-01-01";
                var[:] = wind_time
                
                #data variables
                var = ncout.createVariable("precipitation", float, ("time","lat", "lon"), 
                                           zlib=True, complevel=1, shuffle=True, chunksizes=(1, wind_lat_dimension, wind_lon_dimension));
                var.units = "mm d-1";
                var.long_name = "Pre storm (15 days) mean of ERA5 hourly precipitation on a 0.25X0.25 degree spatial with hourly temporal resolution";
                var[:] = precip_on_wind_grid_prestormref;
                
                ncout.close();   
                
                #Clean up and close datasets and large arrays
                # 1. Close the input NetCDF file (Crucial!)
                if 'precip_nc' in locals() and precip_nc.isopen():
                    precip_nc.close()

                # 2. Delete large 3D Data Arrays (The biggest memory hogs)
                del p_data_raw, p_data
                del precip_on_wind_grid
                del precip_on_wind_grid_prestormref
                
                # 3. Delete calculation intermediates
                del precip_prestormref
                del idx, abs_deltas  # Created inside the loop
                
                # 4. Delete coordinate/time arrays created for precip
                del p_lat, p_lon, p_time, precip_dates
                
                # 5. Delete File Handlers/Paths
                del ncout, precip_nc, processedFilePath, var

                # 6. Force Garbage Collection (Optional but recommended for large loops)
                gc.collect()

                #Print to console that precipitation regridded.
                print("precip regridded for Storm = "+storm)             
                
                                
                #### MAXSS ERA5 pressure data

                #### load data
                pressure_nc = nc.Dataset(path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\MAXSS_HIST_TC_{3}_{1}_{4}_ERA5_SLP.nc".format(region,year,storm,region_id,storm_id)));
                #pressure = pressure_nc.variables['__eo_sp'][:]#Sea level pressure
                pressure_lat = pressure_nc.variables['lat'][:]
                pressure_lon = pressure_nc.variables['lon'][:]
                pressure_time = pressure_nc.variables['time'][:]
                
                # Load pressure variable (__eo_sp)
                # Use np.ma.filled to convert masked values to NaNs immediately
                pressure_data_raw = pressure_nc.variables['__eo_sp'][:] 
                pressure_data = np.ma.filled(pressure_data_raw, fill_value=np.nan)
                
                #### Grid Safety Check (Crucial when skipping regridding)
                # Ensure latitude orientation matches wind (e.g., both descending or both ascending)
                if pressure_lat[0] > pressure_lat[-1] and wind_lat[0] < wind_lat[-1]:
                     pressure_lat = pressure_lat[::-1]
                     pressure_data = pressure_data[:, ::-1, :]
                elif pressure_lat[0] < pressure_lat[-1] and wind_lat[0] > wind_lat[-1]:
                     pressure_lat = pressure_lat[::-1]
                     pressure_data = pressure_data[:, ::-1, :]
                
                # Verify shapes match
                if (wind_lat.shape != pressure_lat.shape) or (wind_lon.shape != pressure_lon.shape):
                    print(f"\n!!! FATAL ERROR: Pressure Grid Dimension Mismatch for {storm} !!!")
                    print(f"  Wind: {wind_lat.shape}, Pressure: {pressure_lat.shape}")
                    sys.exit(1)
                
                # Verify coordinates match (allow for tiny float differences)
                if not np.allclose(wind_lat, pressure_lat, atol=1e-4) or not np.allclose(wind_lon, pressure_lon, atol=1e-4):
                    print(f"\n!!! FATAL ERROR: Pressure Grid Coordinates do not match Wind Grid for {storm} !!!")
                    sys.exit(1)
                    
                #### Vectorized Time Alignment
                # Instead of looping, we align time indices efficiently
                pressure_dates = num2date(pressure_time, pressure_nc.variables['time'].units)
                pressure_on_wind_grid = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=np.float32)

                # Vectorized search for closest time index
                # (This is much faster than the nested loop method)
                for wind_step in range(wind_time_dimension):
                    # Find index of pressure time closest to current wind time
                    # If the time axes are identical, this simply maps 0->0, 1->1, etc.
                    abs_deltas = np.absolute(pressure_dates - wind_dates[wind_step])
                    idx = np.argmin(abs_deltas)
                    pressure_on_wind_grid[wind_step, :, :] = pressure_data[idx, :, :]
                                                
                #### save pressure pre storm output into a netCDF 
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_MAXSS_ERA5_pressure.nc".format(region,year,storm)));

                ncout = Dataset(processedFilePath, 'w');
                    
                #### create dataset and provide dimensions
                ncout.createDimension("lat", wind_lat_dimension);
                ncout.createDimension("lon", wind_lon_dimension);
                ncout.createDimension("time", wind_time_dimension);
                
                #dimension variables
                var = ncout.createVariable("lat", float, ("lat",));
                var.units = "lat (degrees North)";
                var[:] = wind_lat;
                
                var = ncout.createVariable("lon", float, ("lon",));
                var.units = "lon (degrees East)";
                var[:] = wind_lon;
                
                var = ncout.createVariable("time", int, ("time",));
                var.long_name = "Time";
                var.units = "seconds since 1981-01-01";
                var[:] = wind_time
                
                #data variables
                var = ncout.createVariable("sea_level_pressure", float, ("time","lat", "lon"), 
                                           zlib=True, complevel=1, shuffle=True, chunksizes=(1, wind_lat_dimension, wind_lon_dimension));
                var.units = "Pa";
                var.long_name = "ERA5 hourly sea level pressure resampled to a 0.25X0.25 degree spatial and hourly temporal resolution";
                var[:] = pressure_on_wind_grid;
                
                ncout.close();   
                
                #### MAXSS ERA5 pressure data pre storm_reference
                
                #pre storm is first 15 days,so 15 days * hourly resolution
                #Suppress the "Mean of empty slice" for this calculation as some cells legitimatly always land
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    pressure_prestormref = np.nanmean(pressure_on_wind_grid[0:(15*24):1,:,:], axis=0)
                #create empty grid
                pressure_on_wind_grid_prestormref = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                #set all the values equal to the pre storm mean values
                for wind_step in range(0, wind_time_dimension):
                    pressure_on_wind_grid_prestormref[wind_step,:,:]=pressure_prestormref[:,:]
                

                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_MAXSS_ERA5_pressure_pre_storm_reference.nc".format(region,year,storm)));
                ncout = Dataset(processedFilePath, 'w');
                    #### create dataset and provide dimensions
                
                ncout.createDimension("lat", wind_lat_dimension);
                ncout.createDimension("lon", wind_lon_dimension);
                ncout.createDimension("time", wind_time_dimension);
                
                #dimension variables
                var = ncout.createVariable("lat", float, ("lat",));
                var.units = "lat (degrees North)";
                var[:] = wind_lat;
                
                var = ncout.createVariable("lon", float, ("lon",));
                var.units = "lon (degrees East)";
                var[:] = wind_lon;
                
                var = ncout.createVariable("time", int, ("time",));
                var.long_name = "Time";
                var.units = "seconds since 1981-01-01";
                var[:] = wind_time
                
                #data variables
                var = ncout.createVariable("sea_level_pressure", float, ("time","lat", "lon"), 
                                           zlib=True, complevel=1, shuffle=True, chunksizes=(1, wind_lat_dimension, wind_lon_dimension));
                var.units = "Pa";
                var.long_name = "Pre storm (15 days) mean of ERA5 hourly pressure resampled to a 0.25X0.25 degree spatial and hourly temporal resolution";
                var[:] = pressure_on_wind_grid_prestormref;
                
                ncout.close();   
                
                #### Cleanup
                if 'pressure_nc' in locals() and pressure_nc.isopen():
                    pressure_nc.close()
                
                del pressure_data_raw, pressure_data, pressure_on_wind_grid
                del pressure_on_wind_grid_prestormref
                del pressure_lat, pressure_lon, pressure_time, pressure_dates
                
                gc.collect()
                
                print("Pressure processed for Storm = "+storm)

                #### MAXSS Land fraction 
                
                #### create dataset and provide dimensions
                winds_nc = nc.Dataset(path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\MAXSS_HIST_TC_{3}_{1}_{4}_MAXSS_HIST_TC_L4.nc".format(region,year,storm,region_id,storm_id)));
                wind_land_fraction = winds_nc.variables['__eo_land_fraction'][0]
                
                
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_MAXSS_land_fraction.nc".format(region,year,storm)));
                ncout = Dataset(processedFilePath, 'w');

                ncout.createDimension("latitude", wind_lat_dimension);
                ncout.createDimension("longitude", wind_lon_dimension);
                
                #dimension variables
                var = ncout.createVariable("latitude", float, ("latitude",));
                var.units = "lat (degrees North)";
                var[:] = wind_lat;
                
                var = ncout.createVariable("longitude", float, ("longitude",));
                var.units = "lon (degrees East)";
                var[:] = wind_lon;
                
                
                #data variables
                var = ncout.createVariable("land_proportion", float, ("latitude", "longitude"), 
                                           zlib=True, complevel=1, shuffle=True, chunksizes=(wind_lat_dimension, wind_lon_dimension));
                var.units = "1";
                var.long_name = "Fraction of grid cell as land (0-1)";
                var[:] = wind_land_fraction;
                
                ncout.close();   
                
                print("Land fraction for Storm = "+storm)

            
                # Following section extensively updated from Richards code to use
                # Ocean Flux dataset from Dan Ford which is more up to date.
                
                #### Load in monthly pco2 data (Ford et.,al 2025)
                #https://zenodo.org/records/15053676
                
                # Setup the download root
                downloadedRoot = "E:/MAXSS_working_directory/Ford_et_al_GBC_fco2/flux"
                
                # Set up month abbreviations: ['jan', 'feb', 'mar', ..., 'dec']
                month_abbrs = [calendar.month_name[i][:3].lower() for i in range(1, 13)]
                
                # set up download file template
                downloaded_file_template = Template(
                    path.join(downloadedRoot, "${YYYY}", "${MM}", "OceanFluxGHG-month${MM}-${MMM}-${YYYY}-v0.nc"))
                    
                # 1. DEFINE TARGET REGION 
                # wind_lon is your original -180 to 180 array.
                # target_lons_360 is the version we use for interpolation.
                target_lons_360 = (wind_lon + 360) % 360
                target_lats = wind_lat
                
                # Create the 2D meshgrid using the 0-360 longitudes
                region_grid_x, region_grid_y = np.meshgrid(target_lons_360, target_lats)
                
                # 2. PREPARE STORAGE ARRAYS
                # Dimensions derived directly from your wind grid
                subset_lat_dim = len(target_lats)
                subset_lon_dim = len(target_lons_360)
                
                all_pco2_data_region_subset = np.empty((12, subset_lat_dim, subset_lon_dim), dtype=float)
                all_conc_co2_air_data_region_subset = np.empty((12, subset_lat_dim, subset_lon_dim), dtype=float)
                all_reynolds_data_region_subset = np.empty((12, subset_lat_dim, subset_lon_dim), dtype=float)

                # 3. LOOP THROUGH MONTHS
                for imonth in range(0, 12):
                    monthStr = f"{imonth+1:02d}"
                    monthAbbr = month_abbrs[imonth]
                    downloadFilePath = downloaded_file_template.safe_substitute(YYYY=year, MM=monthStr, MMM=monthAbbr)
                    print(f"Processing Month {imonth+1}: {downloadFilePath}")
                
                    try:
                        with nc.Dataset(downloadFilePath, 'r') as pco2_nc:
                            
                            # --- A. DYNAMIC COORDINATE READING ---
                            # Read what is actually in the file (handling common variable names)
                            if 'lon' in pco2_nc.variables:
                                src_lon = pco2_nc.variables['lon'][:]
                                src_lat = pco2_nc.variables['lat'][:]
                            elif 'longitude' in pco2_nc.variables:
                                src_lon = pco2_nc.variables['longitude'][:]
                                src_lat = pco2_nc.variables['latitude'][:]
                            else:
                                raise ValueError("Could not find 'lon' or 'longitude' variable in NetCDF.")

                            # Detect Source Domain (Is it 0-360 or -180 to 180?)
                            source_is_180 = np.min(src_lon) < 0
                            
                            # --- B. ADJUST TARGET GRID TO MATCH SOURCE ---
                            # We adjust our target grid (region_grid_x) to match the file's system
                            # so griddata can find the points.
                            
                            adjusted_region_grid_x = region_grid_x.copy()
                            
                            if source_is_180:
                                # Source is -180 to 180. 
                                # Convert our 0-360 target grid to -180 to 180.
                                # (Points > 180 become negative)
                                adjusted_region_grid_x = np.where(adjusted_region_grid_x > 180, 
                                                                  adjusted_region_grid_x - 360, 
                                                                  adjusted_region_grid_x)
                            else:
                                # Source is 0 to 360.
                                # Ensure our target grid is 0-360
                                adjusted_region_grid_x = adjusted_region_grid_x % 360

                            # --- C. SAFETY CHECK ---
                            if pco2_nc.variables["OBPC"].shape[0] != 1:
                                raise ValueError("CRITICAL ERROR: File contains multiple timesteps.")

                            # Extract variables
                            raw_pco2 = pco2_nc.variables["OBPC"][0, :, :] #this is the raw pco2 loaded in from Dan F (#Uexp-fnn-u) data
                            raw_air = pco2_nc.variables["V_gas"][0, :, :]
                            raw_sst = pco2_nc.variables["FT1_Kelvin_mean"][0, :, :]
                            
                            # --- D. WRAP SOURCE DATA (Seamless Global) ---
                            # Wrap coords to close the global seam
                            src_lon_wrapped = np.concatenate(([src_lon[-1] - 360], src_lon, [src_lon[0] + 360]))
                            
                            # Pad Data
                            pco2_wrapped = np.pad(raw_pco2, ((0,0), (1,1)), mode='wrap')
                            air_wrapped  = np.pad(raw_air,  ((0,0), (1,1)), mode='wrap')
                            sst_wrapped  = np.pad(raw_sst,  ((0,0), (1,1)), mode='wrap')
                            
                            # Flatten Source for griddata
                            src_lon_grid_w, src_lat_grid_w = np.meshgrid(src_lon_wrapped, src_lat)
                            flat_lons = src_lon_grid_w.flatten()
                            flat_lats = src_lat_grid_w.flatten()
                            flat_pco2 = pco2_wrapped.flatten()
                            flat_air  = air_wrapped.flatten()
                            flat_sst  = sst_wrapped.flatten()

                            # --- E. FILTER MASK (Fixes Land Bleeding) ---
                            # Only interpolate using valid ocean points
                            valid_mask = ~np.ma.getmaskarray(flat_pco2) & ~np.isnan(flat_pco2)
                            
                            points_valid = np.column_stack((flat_lons[valid_mask], flat_lats[valid_mask]))
                            values_pco2_valid = flat_pco2[valid_mask]
                            values_air_valid  = flat_air[valid_mask]
                            values_sst_valid  = flat_sst[valid_mask]

                            # --- F. INTERPOLATE ---
                            all_pco2_data_region_subset[imonth, :, :] = griddata(
                                points_valid, 
                                values_pco2_valid, 
                                (adjusted_region_grid_x, region_grid_y), 
                                method='linear'
                            )
                            
                            all_conc_co2_air_data_region_subset[imonth, :, :] = griddata(
                                points_valid, 
                                values_air_valid, 
                                (adjusted_region_grid_x, region_grid_y), 
                                method='linear'
                            )
                            
                            all_reynolds_data_region_subset[imonth, :, :] = griddata(
                                points_valid, 
                                values_sst_valid, 
                                (adjusted_region_grid_x, region_grid_y), 
                                method='linear'
                            )
                            
                            # --- F: POST-INTERPOLATION LAND MASK ---
                            # 1. Build a Tree of your VALID source points (the ocean data)
                            #    (points_valid is already defined in your code as the X/Y of real ocean data)
                            tree = cKDTree(points_valid)
                            
                            # 2. Prepare your Target Grid points for querying
                            #    We need a list of every (Lon, Lat) you just interpolated onto
                            #    (Flattening the 2D meshgrids into a list of points)
                            target_points_flat = np.column_stack((adjusted_region_grid_x.flatten(), region_grid_y.flatten()))
                            
                            # 3. Find the distance to the nearest valid source point for every target point
                            #    dists will be an array of distances in degrees
                            dists, _ = tree.query(target_points_flat)
                            
                            # 4. Reshape the distance array back to the 2D grid shape (lat, lon)
                            dists_grid = dists.reshape(adjusted_region_grid_x.shape)
                            
                            # 5. Apply the Mask
                            #    The Ford dataset is likely 1x1 degree resolution. 
                            #    The diagonal of a 1x1 box is ~1.41 degrees.
                            #    So, if a point is > 1.5 degrees from data, it's likely in a "gap" (land).
                            GAP_THRESHOLD = 1.5  # Degrees. Adjust this if needed (e.g. 2.0 or 1.2)
                            
                            # Create a boolean mask where the distance is too big
                            is_land_gap = dists_grid > GAP_THRESHOLD
                            
                            # Set those points to NaN in your data arrays
                            all_pco2_data_region_subset[imonth, :, :][is_land_gap] = np.nan
                            all_conc_co2_air_data_region_subset[imonth, :, :][is_land_gap] = np.nan
                            all_reynolds_data_region_subset[imonth, :, :][is_land_gap] = np.nan
                            
                    except FileNotFoundError:
                        print(f"Warning: File not found for {monthAbbr}-{year}")
                
                # make pCO2 matrix the same temporal scale as wind data
                
                #create empty pco2 grid the same size as the wind data.
                pco2_on_wind_grid = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                pco2_on_wind_grid[:] = np.nan
                conc_pco2_air_on_wind_grid = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                conc_pco2_air_on_wind_grid[:] = np.nan
                reynolds_co2_on_wind_grid = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                reynolds_co2_on_wind_grid[:] = np.nan
                
                
                #loop through the wind timestamps, extract the month and use that to pick which pco2 data to use.
                for wind_step in range(0, wind_time_dimension):
                
                    a=wind_dates[wind_step].timetuple()
                    #this is month a[1]
                    #need to index at 0 so -1 month
                    month=a[1]-1
                    
                    pco2_on_wind_grid[wind_step,:,:]=all_pco2_data_region_subset[month,:,:]
                    conc_pco2_air_on_wind_grid[wind_step,:,:]=all_conc_co2_air_data_region_subset[month,:,:]
                    reynolds_co2_on_wind_grid[wind_step,:,:]=all_reynolds_data_region_subset[month,:,:]
                
                del all_pco2_data_region_subset
                del all_conc_co2_air_data_region_subset
                del all_reynolds_data_region_subset
                gc.collect()
                
                # 4. save pco2 output into a netCDF 
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_Ford_et_al_pco2.nc".format(region,year,storm)));
                ncout = Dataset(processedFilePath, 'w');
                # create dataset and provide dimensions
                
                ncout.createDimension("lat", wind_lat_dimension);
                ncout.createDimension("lon", wind_lon_dimension);
                ncout.createDimension("time", wind_time_dimension);
                
                #dimension variables
                var = ncout.createVariable("lat", float, ("lat",));
                var.units = "lat (degrees North)";
                var[:] = wind_lat;
                
                var = ncout.createVariable("lon", float, ("lon",));
                var.units = "lon (degrees East)";
                var[:] = wind_lon;
                
                var = ncout.createVariable("time", int, ("time",));
                var.long_name = "Time";
                var.units = "seconds since 1981-01-01";
                var[:] = wind_time
                
                #data variables
                var = ncout.createVariable("pCO2water_mean", float, ("time","lat", "lon"), 
                                           zlib=True, complevel=1, shuffle=True, chunksizes=(1, wind_lat_dimension, wind_lon_dimension));
                var.units = "microatm";
                var.long_name = "Monthly pco2 resampled to hourly on a 0.25X0.25 degree spatial resolution";
                var[:] = pco2_on_wind_grid;
                
                #data variables
                var = ncout.createVariable("V_gas", float, ("time","lat", "lon"), 
                                           zlib=True, complevel=1, shuffle=True, chunksizes=(1, wind_lat_dimension, wind_lon_dimension));
                var.units = "micromol mol-1";
                var.long_name = "Concentration CO2 in dry air (ppm) or dry molecular fraction of CO2 in the atmosphere (umol mol-1)";
                var[:] = conc_pco2_air_on_wind_grid;
                
                #data variables
                var = ncout.createVariable("reynolds_temperature_mean", float, ("time","lat", "lon"), 
                                           zlib=True, complevel=1, shuffle=True, chunksizes=(1, wind_lat_dimension, wind_lon_dimension));
                var.units = "kelvin";
                var.long_name = "Monthly Sea surface skin temperature resampled to hourly on a 0.25X0.25 degree spatial resolution";
                var[:] = reynolds_co2_on_wind_grid;
                
                ncout.close();   
                
                ###############################################################
                ## Create the pre-storm reference file for pco2 data
                
                #calculate pco2 water pre storm ref              
                #pre storm is first 15 days,so 15 days * hourly resolution
                #Suppress the "Mean of empty slice" for this calculation as some cells legitimatly always land
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    pc02_prestormref = np.nanmean(pco2_on_wind_grid[0:(15*24):1,:,:], axis=0)
                #create empty grid
                pco2_on_wind_grid_prestormref = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                #set all the values equal to the pre storm mean values
                for wind_step in range(0, wind_time_dimension):
                    pco2_on_wind_grid_prestormref[wind_step,:,:]=pc02_prestormref[:,:]
                
                #calculate V gas water pre storm ref              
                #pre storm is first 15 days,so 15 days * hourly resolution
                #Suppress the "Mean of empty slice" for this calculation as some cells legitimatly always land
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    v_gas_prestormref = np.nanmean(conc_pco2_air_on_wind_grid[0:(15*24):1,:,:], axis=0)
                #create empty grid
                v_gas_on_wind_grid_prestormref = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                #set all the values equal to the pre storm mean values
                for wind_step in range(0, wind_time_dimension):
                    v_gas_on_wind_grid_prestormref[wind_step,:,:]=v_gas_prestormref[:,:]
                
                ## Create the reynolds temp mean pre-storm ref
                #pre storm is first 15 days,so 15 days * hourly resolution
                #Suppress the "Mean of empty slice" for this calculation as some cells legitimatly always land
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    reynolds_temp_prestormref = np.nanmean(reynolds_co2_on_wind_grid[0:(15*24):1,:,:], axis=0)
                #create empty grid
                reynolds_temp_on_wind_grid_prestormref = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                #set all the values equal to the pre storm mean values
                for wind_step in range(0, wind_time_dimension):
                    reynolds_temp_on_wind_grid_prestormref[wind_step,:,:]=reynolds_temp_prestormref[:,:]
                
                #save pco2 output into a netCDF 
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_Ford_et_al_pco2_pre_storm_reference.nc".format(region,year,storm)));
                ncout = Dataset(processedFilePath, 'w');
                # create dataset and provide dimensions
                
                ncout.createDimension("lat", wind_lat_dimension);
                ncout.createDimension("lon", wind_lon_dimension);
                ncout.createDimension("time", wind_time_dimension);
                
                #dimension variables
                var = ncout.createVariable("lat", float, ("lat",));
                var.units = "lat (degrees North)";
                var[:] = wind_lat;
                
                var = ncout.createVariable("lon", float, ("lon",));
                var.units = "lon (degrees East)";
                var[:] = wind_lon;
                
                var = ncout.createVariable("time", int, ("time",));
                var.long_name = "Time";
                var.units = "seconds since 1981-01-01";
                var[:] = wind_time
                
                #data variables
                var = ncout.createVariable("pCO2water_mean", float, ("time","lat", "lon"), 
                                           zlib=True, complevel=1, shuffle=True, chunksizes=(1, wind_lat_dimension, wind_lon_dimension));
                var.units = "microatm";
                var.long_name = "Pre storm (15 days) mean of monthly pco2 resampled to hourly on a 0.25X0.25 degree spatial resolution";
                var[:] = pco2_on_wind_grid_prestormref;
                
                #data variables
                var = ncout.createVariable("V_gas", float, ("time","lat", "lon"), 
                                           zlib=True, complevel=1, shuffle=True, chunksizes=(1, wind_lat_dimension, wind_lon_dimension));
                var.units = "micromol mol-1";
                var.long_name = "Pre storm (15 days) mean of concentration CO2 in dry air (ppm) or dry molecular fraction of CO2 in the atmosphere (umol mol-1)";
                var[:] = v_gas_on_wind_grid_prestormref;
                
                #data variables
                var = ncout.createVariable("reynolds_temperature_mean", float, ("time","lat", "lon"), 
                                           zlib=True, complevel=1, shuffle=True, chunksizes=(1, wind_lat_dimension, wind_lon_dimension));
                var.units = "kelvin";
                var.long_name = "Pre storm (15 days) mean of monthly Sea surface skin temperature resampled to hourly on a 0.25X0.25 degree spatial resolution";
                var[:] = reynolds_temp_on_wind_grid_prestormref;
                
                ncout.close();   
                
                #delete no longer required variables
                del pco2_on_wind_grid
                del conc_pco2_air_on_wind_grid
                del reynolds_co2_on_wind_grid
                
                del pco2_on_wind_grid_prestormref
                del v_gas_on_wind_grid_prestormref
                del reynolds_temp_on_wind_grid_prestormref
                gc.collect()

                
                print("CO2 regridded for Storm = "+storm)
                
                
                ##NOT SURE IF I NEED TO USE CODE BELOW THIS POINT
                
                
                # #### Verification data - atm pressure data monthly ECMWF
                # downloadedRoot=("verification_data\\air_pressure\\2010");
                # downloadedFileTemplate = Template(path.join(downloadedRoot,"2010${MM}_OCF-PRE-GLO-1M-100-ECMWF.nc"));
                # all_air_pressure = np.empty((12, 180, 360), dtype=float);
                
                # all_air_pressure_region_subset = np.empty((12, wind_lat_dimension, wind_lon_dimension), dtype=float);
                
                
                #     #### resample atmospheric pressure to wind spatial grid
                # for imonth in range(0, 12):
                #     monthStr = format(imonth+1, "02d");
                #     downloadFilePath = downloadedFileTemplate.safe_substitute(MM=monthStr);
                #     air_pressure_nc = nc.Dataset(downloadFilePath, 'r')
                #     all_air_pressure[:,:] = air_pressure_nc.variables["msl_mean"][:,:];
                
                
                #     x=all_air_pressure[imonth,:,:]
                
                #     # to do the interpolation, the grid needs to be as column lists, this does that
                    
                #     #these lines list lat and long of all points as columns
                #     x_coord_range = [i for i in range(0, 360, 1)]
                #     y_coord_range = [i for i in range(0, 180, 1)]
                #     xy_coord = list(itertools.product(x_coord_range, y_coord_range))
                    
                #     #turn 2d data grid of data to column
                #     #air pressure
                #     values = x.flatten(order='F')
                
                #     sample_df = pd.DataFrame()
                
                    
                #     #this is first input into interpolation function
                #     sample_df['X'] = [xy[0] for xy in xy_coord]
                #     sample_df['Y'] = [xy[1] for xy in xy_coord]
                    
                #     #this is second input into interpolation function
                #     sample_df1 = pd.DataFrame()
                #     sample_df1['value'] = values
                
                    
                #     #this is the new grid for the data. global 0.25 degree
                #     grid_x, grid_y = np.mgrid[ 0:360:1440j,0:180:720j]
                
                #     #grid_z0 = griddata(sample_df, sample_df1, (grid_x,grid_y), method='nearest')
                #     grid_z1 = griddata(sample_df, sample_df1, (grid_x,grid_y), method='linear')
                
                
                #     #grid_z2 = griddata(sample_df, sample_df1, (grid_x,grid_y), method='cubic')
                
                #     yyy=grid_z1[:,:,0]
                #     yyy=np.rot90(yyy,1)
                #     yyy=np.flipud(yyy)
                
                
                #     #plt.pcolor(yyy)
                #     lon_new_grid = np.arange(-180, 180, 0.25)
                #     lat_new_grid = np.arange(-90, 90, 0.25)
                
                #     #now work out the subset of the data o extract pco2 for!
                #     max_lat_ind = np.where(lat_new_grid == max_lat)[0][0]
                #     min_lat_ind = np.where(lat_new_grid == min_lat)[0][0]-1
                
                #     max_lon_ind = np.where(lon_new_grid == max_lon)[0][0]
                #     min_lon_ind = np.where(lon_new_grid == min_lon)[0][0]-1
                    
                #     yyy_subset=yyy[min_lat_ind:max_lat_ind,min_lon_ind:max_lon_ind]
                
                #     # plt.pcolor(yyy_subset)
                #     # plt.pcolor(all_air_pressure_region_subset[1,:,:])
                
                #     all_air_pressure_region_subset[imonth,:,:] = yyy_subset;
                
                #     #### make pCO2 matrix the same temporal scale as wind data
                
                # #create empty pco2 grid the same size as the wind data.
                # air_pressure_on_wind_grid = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                
                
                #     #### loop through the wind timestamps, 
                #     #extract the month and use that to pick which pressure data to use.
                # for wind_step in range(0, wind_time_dimension):
                
                #     a=wind_dates[wind_step].timetuple()
                #     #this is month a[1]
                #     #indexing starts at 0 so minus 1
                #     month=a[1]-1
                    
                #     air_pressure_on_wind_grid[wind_step,:,:]=all_air_pressure_region_subset[month,:,:]
                
                #     #### save air pressure output into a netCDF 
                # processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_verification_data_ECMWF_air_pressure.nc".format(region,year,storm)));
                # ncout = Dataset(processedFilePath, 'w');
                #     #### create dataset and provide dimensions
                
                # ncout.createDimension("lat", wind_lat_dimension);
                # ncout.createDimension("lon", wind_lon_dimension);
                # ncout.createDimension("time", wind_time_dimension);
                
                # #dimension variables
                # var = ncout.createVariable("lat", float, ("lat",));
                # var.units = "lat (degrees North)";
                # var[:] = wind_lat;
                
                # var = ncout.createVariable("lon", float, ("lon",));
                # var.units = "lon (degrees East)";
                # var[:] = wind_lon;
                
                # var = ncout.createVariable("time", int, ("time",));
                # var.long_name = "Time";
                # var.units = "seconds since 1981-01-01";
                # var[:] = wind_time
                
                # #data variables
                # var = ncout.createVariable("sea_level_pressure", float, ("time","lat", "lon"));
                # var.units = "pascals";
                # var.long_name = "Mean monthly air pressure resampled to hourly on a 0.25X0.25 degree spatial resolution";
                # var[:] = air_pressure_on_wind_grid;
                
                # ncout.close();   
                # print("Air pressure regridded for Storm = "+storm)

                
 
                
                # #### world seas mask
                # downloadedRoot=("verification_data\\");
                # downloadFilePath = path.join(downloadedRoot,"World_Seas-final-complete_IGA.nc");
                # world_seas = np.empty((180, 360), dtype=float);
                # world_seas_nc = nc.Dataset(downloadFilePath, 'r')
                # world_seas[:,:] = world_seas_nc.variables["sea-mask"][:,:];
                # x=world_seas[:,:]

                
                # #these lines list lat and long of all points as columns
                # x_coord_range = [i for i in range(0, 360, 1)]
                # y_coord_range = [i for i in range(0, 180, 1)]
                # xy_coord = list(itertools.product(x_coord_range, y_coord_range))
                
                # #turn 2d data grid of data to column
                # values = x.flatten(order='F')
                 
                # sample_df = pd.DataFrame()

                # #this is first input into interpolation function
                # sample_df['X'] = [xy[0] for xy in xy_coord]
                # sample_df['Y'] = [xy[1] for xy in xy_coord]
                
                # #this is second input into interpolation function
                # sample_df1 = pd.DataFrame()
                # sample_df1['value'] = values
                
                # values[values == -999] = 'nan' # or use np.nan
        
                # #this is the new grid for the data. global 0.25 degree
                # grid_x, grid_y = np.mgrid[ 0:360:1440j,0:180:720j]
            
                # #grid_z0 = griddata(sample_df, sample_df1, (grid_x,grid_y), method='nearest')
                # #grid_z1 = griddata(sample_df, sample_df1, (grid_x,grid_y), method='linear')
                # grid_z1 = griddata(sample_df, sample_df1, (grid_x,grid_y), method='nearest')

                # #grid_z2 = griddata(sample_df, sample_df1, (grid_x,grid_y), method='cubic')
            
                # yyy=grid_z1[:,:,0]
                # yyy=np.rot90(yyy,1)
                # #yyy=np.flipud(yyy)
            
                # #plt.pcolor(yyy)
                # lon_new_grid = np.arange(-180, 180, 0.25)
                # lat_new_grid = np.arange(-90, 90, 0.25)
            
                # #now work out the subset of the data to extract for!
                # max_lat_ind = np.where(lat_new_grid == max_lat)[0][0]
                # min_lat_ind = np.where(lat_new_grid == min_lat)[0][0]-1
            
                # max_lon_ind = np.where(lon_new_grid == max_lon)[0][0]
                # min_lon_ind = np.where(lon_new_grid == min_lon)[0][0]-1
                
                # yyy_subset=yyy[min_lat_ind:max_lat_ind,min_lon_ind:max_lon_ind]

            

                # #### save mask output into a netCDF 
                # processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_world_ocean_mask.nc".format(region,year,storm)));
                # ncout = Dataset(processedFilePath, 'w');
                # # create dataset and provide dimensions
                
                # ncout.createDimension("lat", wind_lat_dimension);
                # ncout.createDimension("lon", wind_lon_dimension);
                
                # #dimension variables
                # var = ncout.createVariable("lat", float, ("lat",));
                # var.units = "lat (degrees North)";
                # var[:] = wind_lat;
                
                # var = ncout.createVariable("lon", float, ("lon",));
                # var.units = "lon (degrees East)";
                # var[:] = wind_lon;
                

                # #data variables
                # var = ncout.createVariable("seamask", float, ("lat", "lon"));
                # var.units = "unitless";
                # var.long_name = "World seas mask nearest interpolation on to  a 0.25X0.25 degree spatial resolution";
                # var[:] = yyy_subset;
                
                # ncout.close();   
                
                # print("mask regridded for Storm = "+storm)

                
                
                

# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 12:12:13 2026

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
import datetime

import argparse #DJF 09/05/2026: So we can add code to allow this to run from a bash script


def process_slice(valData, errData,wind_lat_dimension, wind_lon_dimension, min_lon, max_lon, min_lat, max_lat, iCoordMeshes, outputRes=0.25):

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


# if __name__ == "__main__": # DJF: 09/05/2026: Turning the _main_ into a function
def MAXSS_resample_main(MAXSS_working_directory = "E:/MAXSS_working_directory", downloadedRoot = "E:/MAXSS_working_directory/Ford_et_al_GBC_fco2/flux", specified_storms = [], MAXSS_regions = ["north-atlantic"], specified_years =[]): # Values after the equals are the default values that would be called if you used function as: MAXSS_resample_main()


    """
    Function to run the MAXSS resampling script and allowing it to be callable.

    Inputs:
    MAXSS_working_directory = Sets the working directory for the script
    downloadedRoot = Sets the directory for the UExP-FNN-U downloaded 'flux' folder
    """
    #### Get the path of the root directory where the data are stored.
    #This will be user specific and can be changed depening on where data is stored
    # MAXSS_working_directory = "E:/MAXSS_working_directory"; #DJF 09/05/2026: this is now brought in wiht the function arguements

    #note to use the same file structure used by the project r.g. #maxss/storm-atlas/ibtracts/region/year/storm
    os.chdir(MAXSS_working_directory);
    print("Working directory is now:", os.getcwd());

    #### Loop through the regions in MAXSS storm dataset
    for region in MAXSS_regions:

        #define the directory for the region
        region_directory = os.path.join(MAXSS_working_directory, "maxss", "storm-atlas", "tropical", "ibtracs", region)
        #look for the year subfolders

        #get a list of the years
        year_list = []
        for entry_name in os.listdir(region_directory):
            entry_path = path.join(region_directory, entry_name)
            if path.isdir(entry_path):
                year_list.append(entry_name)
        #get a list of the paths for each year folder
        year_directory_list = [os.path.join(region_directory, year) for year in year_list]

        #define to loop through years
        MAXSS_years=year_list

        #### Loop through the years and their directory paths simultaneously
        # zip() pairs up year_list ["2010", "2011"] with year_directory_list ["/path/2010", "/path/2011"]
        for year, current_year_dir in zip(MAXSS_years, year_directory_list):
            
            # 2. Add the filtering block for years right here:
            if len(specified_years) > 0:
                # If the current year doesn't match any fragment in specified_years, skip it!
                if not any(str(yr_fragment) in year for yr_fragment in specified_years):
                    print(f"Skipping year (not in specified list): {year}")
                    continue

            print(f"Entering Year Directory: {year}")
                
            # Make a list of the storm folder names
            storm_list = []
            for entry_name in os.listdir(current_year_dir):
                entry_path = path.join(current_year_dir, entry_name)
                if path.isdir(entry_path):
                    storm_list.append(entry_name)

            # Build storm directory list
            storm_directory_list = [os.path.join(current_year_dir, storm) for storm in storm_list]

            #define to loop through years
            MAXSS_storms=storm_list

            #### Loop through the storms and their directory paths simultaneously using zip()
            for storm, current_storm_dir in zip(MAXSS_storms, storm_directory_list):

                # If specified_storms has entries, check if ANY of your fragments match the storm name.
                if len(specified_storms) > 0:
                    # This checks if NONE of the fragments are inside the storm string
                    if not any(fragment in storm for fragment in specified_storms):
                        print(f"Skipping storm (no fragment match): {storm}")
                        continue  # No manual counters to update anymore! Safe and clean.

                print(f"--> Processing storm: {storm}")

                # Directory for storm being processed
                storm_dir = current_storm_dir
                
                #directory for storm being processes relative to current working directory
                storm_dir_relative = os.path.join("maxss", "storm-atlas", "tropical", "ibtracs", region, year, storm)
                
                output_location = os.path.join(MAXSS_working_directory, "output", "maxss", "storm-atlas", "tropical", "ibtracs", region, year, storm)
                
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
                winds_nc = nc.Dataset(os.path.join(MAXSS_working_directory, "maxss", "storm-atlas", "tropical", "ibtracs", region, year, storm, 
                                                   f"MAXSS_HIST_TC_{region_id}_{year}_{storm_id}_MAXSS_HIST_TC_L4.nc"))

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
                processedFilePath = path.join(MAXSS_working_directory, "maxss", "storm-atlas", "tropical", "ibtracs", region, year, storm, "Resampled_for_fluxengine_storm_timings_with_masks.nc")
                
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
                processedFilePath = path.join(MAXSS_working_directory, "maxss", "storm-atlas", "tropical", "ibtracs", region, year, storm, "Resampled_for_fluxengine_MAXSS_L4_windspeed.nc")

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

                ncout.close()

                gc.collect()

                # 8. Close wind data netcdf
                winds_nc.close()

                print("Wind regridded for Storm = "+storm)


                ## MAXSS ESACCI SST data ##
                # load data
                sst_nc = nc.Dataset(path.join(MAXSS_working_directory, "maxss", "storm-atlas", "tropical", "ibtracs", region, year, storm, f"MAXSS_HIST_TC_{region_id}_{year}_{storm_id}_OSTIA-C3S-L4-GLOB-v2.1.nc"))

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
                    #newVals, newCountCount, newValsErr = process_slice(sst_time_slice, sst_uncertainty_slice);
                    
                    newVals, newCountCount, newValsErr = process_slice(sst_time_slice, sst_uncertainty_slice, 
                        wind_lat_dimension, wind_lon_dimension, min_lon, max_lon, min_lat, max_lat, iCoordMeshes)
                    

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
                processedFilePath = path.join(MAXSS_working_directory, "maxss", "storm-atlas", "tropical", "ibtracs", region, year, storm, "Resampled_for_fluxengine_MAXSS_ESACCI_SST.nc")
                
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
                                           zlib=True, complevel=1, shuffle=True, fill_value=sst_fill_value,
                                           chunksizes=(1, wind_lat_dimension, wind_lon_dimension));
                var.units = "Degrees Kelvin";
                var.long_name = "Daily ESACCI sea surface temperature resampled to a 0.25X0.25 degree spatial and hourly temporal resolutionn";
                var[:] = sst_on_wind_grid;

                ncout.close();

                # . close sst netcdf
                sst_nc.close()


                ## MAXSS ESACCI SSS data ##
                # 1. load data
                sss_nc = nc.Dataset(path.join(MAXSS_working_directory, "maxss", "storm-atlas", "tropical", "ibtracs", region, year, storm, f"MAXSS_HIST_TC_{region_id}_{year}_{storm_id}_ESACCI-SEASURFACESALINITY-L4-SSS-MERGED_OI_7DAY_RUNNINGMEAN_DAILY_25km.nc"))
                                    
                #sss = sss_nc.variables['__eo_sss'][:]
                sss_lat = sss_nc.variables['lat'][:]
                sss_lon = sss_nc.variables['lon'][:]
                sss_time = sss_nc.variables['time'][:]

                # 2. Extract sss fill value
                sss_fill_value = sss_nc.variables['__eo_sss']._FillValue

                # 3. resample SSS to wind grid
                iCoordMeshes = np.full((wind_lat_dimension, wind_lon_dimension), None, dtype=object)

                CCILats = sss_nc.variables['lat'][:]
                CCILons = sss_nc.variables['lon'][:]
                
                # Build explicit grid edges using linspace (eliminates np.arange floating-point drift)
                lat_edges = np.linspace(min_lat, max_lat, wind_lat_dimension + 1)
                lon_edges = np.linspace(min_lon, max_lon, wind_lon_dimension + 1)
                
                # Digitize map input coordinates into integer bin indices deterministically
                lat_inds = np.digitize(CCILats, lat_edges) - 1
                lon_inds = np.digitize(CCILons, lon_edges) - 1
                
                for ilat in range(wind_lat_dimension):
                    wlat = np.where(lat_inds == ilat)[0]
                    if len(wlat) == 0:
                        continue
                    for ilon in range(wind_lon_dimension):
                        wlon = np.where(lon_inds == ilon)[0]
                        if len(wlon) > 0:
                            iCoordMeshes[ilat, ilon] = np.meshgrid(wlat, wlon)

                # 4. PRE-COMPUTE COORDINATES FOR CKDTREE
                grd_lons, grd_lats = np.meshgrid(wind_lon, wind_lat)
                target_coords = np.column_stack((grd_lons.ravel(), grd_lats.ravel()))

                # 5. loop through timesteps
                timesteps_sss = len(sss_time)

                # Store data for each day of the month
                sss_regrid_Vals = np.empty((timesteps_sss, wind_lat_dimension, wind_lon_dimension), dtype=float)

                for sss_timesteps in range(0, timesteps_sss):
                    sss_time_slice = sss_nc.variables['__eo_sss'][sss_timesteps, :, :]
                    sss_uncertainty_slice = sss_nc.variables['__eo_sss_random_error'][sss_timesteps, :, :]
                    
                    # Initial binning
                    newVals, newCountCount, newValsErr = process_slice(
                        sss_time_slice, sss_uncertainty_slice, 
                        wind_lat_dimension, wind_lon_dimension, 
                        min_lon, max_lon, min_lat, max_lat, iCoordMeshes
                    )

                    # --- REPLACED GRIDDATA WITH CKDTREE ---
                    valid_mask = ~np.isnan(newVals)
                    valid_lats = grd_lats[valid_mask]
                    valid_lons = grd_lons[valid_mask]
                    valid_values = newVals[valid_mask]

                    if len(valid_values) > 0:
                        source_coords = np.column_stack((valid_lons, valid_lats))
                        
                        # Build spatial tree of valid binned data
                        tree = cKDTree(source_coords)
                        
                        # Query nearest neighbors up to 0.35 degrees (~35km search radius)
                        # Anything further is kept as NaN to prevent coastal/cross-track bleeding
                        max_distance = 0.35  
                        distances, indices = tree.query(target_coords, distance_upper_bound=max_distance)

                        nan_mask = indices == len(source_coords)
                        indices[nan_mask] = 0  # Avoid IndexErrors
                        
                        resampled_flat = valid_values[indices]
                        resampled_flat[nan_mask] = np.nan

                        SSS_nearest = resampled_flat.reshape(wind_lat_dimension, wind_lon_dimension)
                    else:
                        SSS_nearest = np.full((wind_lat_dimension, wind_lon_dimension), np.nan)

                    # Apply SST land mask (your existing mask logic)
                    where_sst_nan = np.isnan(sst_regrid_Vals[0, :, :])
                    SSS_nearest[where_sst_nan] = np.nan

                    sss_regrid_Vals[sss_timesteps, :, :] = SSS_nearest

                # 6. Get the data of each timestep in wind data
                sss_on_wind_grid = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float)
                sss_time = sss_nc.variables['time'][:]
                sss_dates = num2date(sss_time, sss_nc.variables['time'].units)

                # 7. loop through the wind timestamps, extract the dates and find the closest SST data.
                for wind_step in range(0, wind_time_dimension):
                    abs_deltas_from_target_date = np.absolute(sss_dates - wind_dates[wind_step])
                    index_of_min_delta_from_target_date = np.argmin(abs_deltas_from_target_date)

                    sss_on_wind_grid[wind_step, :, :] = sss_regrid_Vals[index_of_min_delta_from_target_date, :, :]

                # 8. Apply the 3D analysis mask
                sss_on_wind_grid[analysis_period_mask_3d == 0] = np.nan

                # 9. Convert all NaNs (land, mask ect.) to the SSS fill value
                sss_on_wind_grid = np.nan_to_num(sss_on_wind_grid, nan=sss_fill_value).astype('float32')

                # 10. Save SSS output into a netCDF
                processedFilePath = path.join(MAXSS_working_directory, "maxss", "storm-atlas", "tropical", "ibtracs", region, year, storm, "Resampled_for_fluxengine_MAXSS_ESACCI_SSS.nc")
   
                ncout = Dataset(processedFilePath, 'w')

                # create dataset and provide dimensions
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
                var = ncout.createVariable("sss", float, ("time", "lat", "lon"),
                                           zlib=True, complevel=1, shuffle=True, fill_value=sss_fill_value, chunksizes=(1, wind_lat_dimension, wind_lon_dimension))
                var.units = "PSU"
                var.long_name = "Daily ESACCI sea surface salinity resampled to a 0.25X0.25 degree spatial and hourly temporal resolution"
                var[:] = sss_on_wind_grid

                ncout.close()

                # 11. close salinity netcdf
                sss_nc.close()

                print("SSS regridded for Storm = " + storm)
                print(storm + " RESAMPLING COMPLETED")

# if __name__ == "__main__":
#     description = """Runs the MAXSS resampling for each storm""";
#     parser = argparse.ArgumentParser(description=description);
#     parser.add_argument("--MAXSS_working_directory", type=str, default="E:/MAXSS_working_directory",
#                         help="Path to the required working directory");
#     parser.add_argument("--downloadedRoot", type=str, default="E:/MAXSS_working_directory/Ford_et_al_GBC_fco2/flux",
#                         help="Path to the data directory containing the UExP-FNN-U flux folder");
    
#     clArgs = parser.parse_args();
#     MAXSS_resample_main(MAXSS_working_directory = clArgs.MAXSS_working_directory, downloadedRoot = clArgs.downloadedRoot)

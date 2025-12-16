# -*- coding: utf-8 -*-
"""

Script for looping through the MAXSS storm-Atlas and resampling data for fluxengine runs

@Original author: Richard Sims

This script has been updated in 2026 by Daniel Wilson

"""

#### all functions and scripts imported here

#Install required packages
import os
from os import path, makedirs;
from glob import glob
from pathlib import Path
import shutil;
import netCDF4 as nc
from netCDF4 import date2num, num2date, Dataset
import argparse;
from string import Template;
import numpy as np;
import urllib.request as request;
from contextlib import closing;
import ssl;
import urllib;
import calendar;
from datetime import datetime, timedelta;
import tarfile;
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import itertools
import pandas as pd
from scipy.interpolate import griddata

def process_slice(valData, errData, outputRes=0.25):
    
    #Handle masked values explicitly to avoid masked element to nan error
    valData = valData.filled(fill_value=np.nan)
    errData = errData.filled(fill_value=np.nan)
    
    newGrid = np.full((wind_lat_dimension, wind_lon_dimension), np.nan, dtype=float);
    newGridCount = np.zeros((wind_lat_dimension, wind_lon_dimension), dtype=float);
    newGridErr = np.full((wind_lat_dimension, wind_lon_dimension), np.nan, dtype=float);
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
            entry_path = os.path.join(region_directory, entry_name)
            if os.path.isdir(entry_path):
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
                entry_path = os.path.join(year_directory_list[year_counter], entry_name)
                if os.path.isdir(entry_path):
                    storm_list.append(entry_name)
            #get a list of the paths for each year folder
            storm_directory_list=glob(year_directory_list[year_counter]+"/*/", recursive = True)
            year_counter=year_counter+1
            
            #define to loop through years
            MAXSS_storms=storm_list
            
            storm_counter=0
            #### Loop through the storms for each year in the MAXSS storm dataset 
            for storm in MAXSS_storms:

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
                    
                #### calculate wind speed from east and west components
                wind_speed = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=np.float32);
                wind_speed=pow((pow(wind_eastward,2)+pow(wind_northward,2)),0.5)
                
                Wind_moment2=pow(wind_speed,2)
                Wind_moment3=pow(wind_speed,3)
                Wind_moment3point7=pow(wind_speed,3.742)

                #### stop wind fields clogging up memory
                del wind_eastward, wind_northward ,winds_nc
                
                #### save wind speed to netCDF
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_MAXSS_L4_windspeed.nc".format(region,year,storm)));

                ncout = Dataset(processedFilePath, 'w');
                
                #### provide dimensions
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
                var = ncout.createVariable("windspeed", float, ("time","lat", "lon"));
                var.units = "ms-1";
                var.long_name = "Hourly wind speed from MAXSS on a 0.25X0.25 degree gridspatial resolution";
                var[:] = wind_speed;
                
                #data variables
                var = ncout.createVariable("second_moment_wind", float, ("time","lat", "lon"));
                var.units = "ms-1 - CHECK";
                var.long_name = "Second moment of wind speed from MAXSS on a 0.25X0.25 degree gridspatial resolution";
                var[:] = Wind_moment2;
                
                #data variables
                var = ncout.createVariable("third_moment_wind", float, ("time","lat", "lon"));
                var.units = "ms-1 - CHECK ";
                var.long_name = "Third moment of wind speed from MAXSS on a 0.25X0.25 degree gridspatial resolution";
                var[:] = Wind_moment3;
                
                var = ncout.createVariable("thirdseven_moment_wind", float, ("time","lat", "lon"));
                var.units = "ms-1 - CHECK";
                var.long_name = "3.7 moment of wind speed from MAXSS on a 0.25X0.25 degree gridspatial resolution";
                var[:] = Wind_moment3point7;
                ncout.close();  
                
                
                #### MAXSS L4 Wind data pre storm reference
                
                #pre storm is first 15 days,so 15 days * hourly resolution
                wind_prestormref=np.nanmean(wind_speed[0:(15*24):1,:,:],axis =(0))
                wind_prestormref_second_moment=np.nanmean(Wind_moment2[0:(15*24):1,:,:],axis =(0))
                wind_prestormref_third_moment=np.nanmean(Wind_moment3[0:(15*24):1,:,:],axis =(0))
                wind_prestormref_thirdseven_moment=np.nanmean(Wind_moment3point7[0:(15*24):1,:,:],axis =(0))

                #create empty grid
                wind_speed_prestormref = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                second_moment_prestormref = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                third_moment_prestormref = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                thirdseven_moment_prestormref = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);

                
                #set all the values equal to the pre storm mean values
                for wind_step in range(0, wind_time_dimension):
                    wind_speed_prestormref[wind_step,:,:]=wind_prestormref[:,:]
                    second_moment_prestormref[wind_step,:,:]=wind_prestormref_second_moment[:,:]
                    third_moment_prestormref[wind_step,:,:]=wind_prestormref_third_moment[:,:]
                    thirdseven_moment_prestormref[wind_step,:,:]=wind_prestormref_thirdseven_moment[:,:]


                #### save wind pre storm output into a netCDF 
                
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_MAXSS_L4_windspeed_pre_storm_reference.nc".format(region,year,storm)));

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
                var = ncout.createVariable("windspeed", float, ("time","lat", "lon"));
                var.units = "ms-1";
                var.long_name = "Pre storm (15 days) mean of hourly wind speed from MAXSS on a 0.25X0.25 degree spatial at a hourly temporal resolution";
                var[:] = wind_speed_prestormref;
                
                #data variables
                var = ncout.createVariable("second_moment_wind", float, ("time","lat", "lon"));
                var.units = "ms-1 - CHECK";
                var.long_name = "Second moment of wind speed from MAXSS on a 0.25X0.25 degree gridspatial resolution";
                var[:] = second_moment_prestormref;
                
                #data variables
                var = ncout.createVariable("third_moment_wind", float, ("time","lat", "lon"));
                var.units = "ms-1 - CHECK ";
                var.long_name = "Third moment of wind speed from MAXSS on a 0.25X0.25 degree gridspatial resolution";
                var[:] = third_moment_prestormref;
                
                var = ncout.createVariable("thirdseven_moment_wind", float, ("time","lat", "lon"));
                var.units = "ms-1 - CHECK";
                var.long_name = "3.7 moment of wind speed from MAXSS on a 0.25X0.25 degree gridspatial resolution";
                var[:] = thirdseven_moment_prestormref;
                ncout.close();   
                
                del wind_speed,Wind_moment2,Wind_moment3,Wind_moment3point7
                del third_moment_prestormref,thirdseven_moment_prestormref,second_moment_prestormref,wind_speed_prestormref
                del wind_prestormref_second_moment,wind_prestormref_third_moment,wind_prestormref_thirdseven_moment
                del wind_prestormref
                print("Wind regridded for Storm = "+storm)


                
                #### MAXSS ESACCI SST data
                #### load data
                sst_nc = nc.Dataset(path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\MAXSS_HIST_TC_{3}_{1}_{4}_OSTIA-C3S-L4-GLOB-v2.1.nc".format(region,year,storm,region_id,storm_id)));
                #sst = sst_nc.variables['__eo_analysed_sst'][:]
                sst_lat = sst_nc.variables['lat'][:]
                sst_lon = sst_nc.variables['lon'][:]
                sst_time = sst_nc.variables['time'][:]
                
                #### resample SST to wind grid
                iCoordMeshes = None; #Initially this is None but will be calculated exactly once. It's a long calculation so doesn't want to be repeated.
                outputRes = 0.25;
                #### Calculate binning information
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
                
                #### Store data for each day of the month
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
                
                    #### get the data of each timestep in wind data
                sst_time = sst_nc.variables['time'][:]
                sst_dates = num2date(sst_time, sst_nc.variables['time'].units)
                
                    #### loop through the wind timestamps, extract the month and use that to pick which SST data to use.
                for wind_step in range(0, wind_time_dimension):
                
                    #now find the index of the CLOSEST value and use that
                    abs_deltas_from_target_date = np.absolute(sst_dates - wind_dates[wind_step])
                    index_of_min_delta_from_target_date = np.argmin(abs_deltas_from_target_date)
                    #closest_date = sst_dates[index_of_min_delta_from_target_date]

                    sst_on_wind_grid[wind_step,:,:]=sst_regrid_Vals[index_of_min_delta_from_target_date,:,:]
                    
                #### save SST output into a netCDF 
                
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_MAXSS_ESACCI_SST.nc".format(region,year,storm)));

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
                var = ncout.createVariable("sst", float, ("time","lat", "lon"));
                var.units = "Degrees Kelvin";
                var.long_name = "Daily ESACCI sea surface temperature resampled to a 0.25X0.25 degree spatial and hourly temporal resolutionn";
                var[:] = sst_on_wind_grid;
                
                ncout.close();   
                
                
                #### MAXSS ESACCI SST data pre_storm_reference
            
                #example slicing for numpy 3d matrixes
                # yy= np.random.randint(0, 100, size=(4, 5, 6))
                # yyy=np.nanmean(yy[0:3:1,:,:],axis =(0)) 
                #temp=sst_on_wind_grid[0:(15*24):1,:,:]
                
                #pre storm is first 15 days,so 15 days * hourly resolution
                SST_prestormref=np.nanmean(sst_on_wind_grid[0:(15*24):1,:,:],axis =(0))
                #create empty grid
                sst_on_wind_grid_prestormref = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                #set all the values equal to the pre storm mean values
                for wind_step in range(0, wind_time_dimension):
                    sst_on_wind_grid_prestormref[wind_step,:,:]=SST_prestormref[:,:]
                
                
                    
                #### save SST pre storm output into a netCDF 
                
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_MAXSS_ESACCI_SST_pre_storm_reference.nc".format(region,year,storm)));

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
                var = ncout.createVariable("sst", float, ("time","lat", "lon"));
                var.units = "Degrees Kelvin";
                var.long_name = "Pre storm (15 days) mean of daily ESACCI sea surface temperature resampled to a 0.25X0.25 degree spatial and hourly temporal resolution";
                var[:] = sst_on_wind_grid_prestormref;
                
                ncout.close();   
                

                #### delete all the variables used during import and saving of SST      
                del SST_prestormref,sst_on_wind_grid_prestormref,ncout, sst_dates, sst_lat,sst_lon,sst_nc,sst_on_wind_grid
                del newVals, sst_regrid_ValsErr,sst_regrid_ValsCounts,newCountCount,newValsErr
                del sst_time,sst_time_slice,sst_timesteps,sst_uncertainty_slice,
                del processedFilePath,timesteps_sst,abs_deltas_from_target_date,iCoordMeshes,index_of_min_delta_from_target_date
                print("SST regridded for Storm = "+storm)
          
                    


                            
                #### MAXSS ESACCI SSS data
                    #### load data                
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
                var = ncout.createVariable("sss", float, ("time","lat", "lon"));
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
                var = ncout.createVariable("sss", float, ("time","lat", "lon"));
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
                
                #### load data
                precip_nc = nc.Dataset(path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\MAXSS_HIST_TC_{3}_{1}_{4}_ERA5_Precipitation.nc".format(region,year,storm,region_id,storm_id)));
                #precip = precip_nc.variables['__eo_tp'][:]#total precipiation
                precip_lat = precip_nc.variables['lat'][:]
                precip_lon = precip_nc.variables['lon'][:]
                precip_time = precip_nc.variables['time'][:]
                
                #### resample precip to wind grid
                iCoordMeshes = None; #Initially this is None but will be calculated exactly once. It's a long calculation so doesn't want to be repeated.
                outputRes = 0.25;
                
                #### Calculate binning information
                #Only do this once because it's computationally expensive but the same for all time steps
                if iCoordMeshes is None: 
                    CCILats = precip_nc.variables['lat'][:]
                    CCILons = precip_nc.variables['lon'][:]
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
                timesteps_precip=len(precip_time)
                del ilat, ilon, CCILats,CCILons, lat , lon
                
                #### Store data for each day of the month
                precip_regrid_Vals = np.empty((timesteps_precip, wind_lat_dimension,wind_lon_dimension), dtype=float);
                # precip_regrid_ValsErr = np.empty((timesteps_precip, wind_lat_dimension,wind_lon_dimension), dtype=float);
                # precip_regrid_ValsCounts = np.empty((timesteps_precip, wind_lat_dimension,wind_lon_dimension), dtype=float);
                               
                for precip_timesteps in range(0, timesteps_precip): 
                    #print(precip_timesteps)                               
                    precip_time_slice=precip_nc.variables['__eo_tp'][precip_timesteps,:,:]   
                    # no uncertainty data but is needed for function so just use the precipitation instead 
                    # and dont use uncertainty.                       
                    precip_uncertainty_slice=precip_nc.variables['__eo_tp'][precip_timesteps,:,:]                         
                    newVals, newCountCount, newValsErr = process_slice(precip_time_slice,precip_uncertainty_slice);
                
                    precip_regrid_Vals[precip_timesteps,:,:] = newVals;
                    # precip_regrid_ValsErr[precip_timesteps,:,:] = newValsErr;
                    # precip_regrid_ValsCounts[precip_timesteps,:,:] = newCountCount;
                
                #plt.pcolor(precip_regrid_Vals[1,:,:])
                             
                precip_on_wind_grid = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                
                #### get the data of each timestep in wind data
                precip_time = precip_nc.variables['time'][:]
                precip_dates = num2date(precip_time, precip_nc.variables['time'].units)
                
                #### loop through the wind timestamps, extract the month and use that to pick which precip data to use.
                for wind_step in range(0, wind_time_dimension):
                
                    #now find the index of the CLOSEST value and use that
                    abs_deltas_from_target_date = np.absolute(precip_dates - wind_dates[wind_step])
                    index_of_min_delta_from_target_date = np.argmin(abs_deltas_from_target_date)
                    #closest_date = precip_dates[index_of_min_delta_from_target_date]

                    precip_on_wind_grid[wind_step,:,:]=precip_regrid_Vals[index_of_min_delta_from_target_date,:,:]
                                    
                                
                # flux engine expects rain in mm d-1 whereas ERA5 are m every hour
                
                unit_conversion_factor_rain= (24/1)*100
                precip_on_wind_grid=precip_on_wind_grid*unit_conversion_factor_rain
                          
                #### save precipitation pre storm output into a netCDF 
                
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
                var = ncout.createVariable("precipitation", float, ("time","lat", "lon"));
                var.units = "m";
                var.long_name = "ERA5 hourly precipitation resampled to a 0.25X0.25 degree spatial and hourly temporal resolution";
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
                var = ncout.createVariable("precipitation", float, ("time","lat", "lon"));
                var.units = "mm d-1";
                var.long_name = "Pre storm (15 days) mean of ERA5 hourly precipitation resampled to a 0.25X0.25 degree spatial and hourly temporal resolution";
                var[:] = precip_on_wind_grid_prestormref;
                
                ncout.close();   
                
                
                    #### delete all the variables used during import and saving of precip      
                del newVals, newValsErr ,precip_on_wind_grid,ncout,precip_dates, precip_lat,precip_lon,precip_nc
                del precip_time,precip_time_slice,precip_timesteps,precip_uncertainty_slice,
                del processedFilePath,timesteps_precip,abs_deltas_from_target_date,iCoordMeshes,index_of_min_delta_from_target_date
                del  precip_regrid_Vals, newCountCount
                print("precip regridded for Storm = "+storm)                
                
                                
                #### MAXSS ERA5 pressure data

                    #### load data
                pressure_nc = nc.Dataset(path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\MAXSS_HIST_TC_{3}_{1}_{4}_ERA5_SLP.nc".format(region,year,storm,region_id,storm_id)));
                #pressure = pressure_nc.variables['__eo_sp'][:]#Sea level pressure
                pressure_lat = pressure_nc.variables['lat'][:]
                pressure_lon = pressure_nc.variables['lon'][:]
                pressure_time = pressure_nc.variables['time'][:]
                
                    #### resample pressure to wind grid
                iCoordMeshes = None; #Initially this is None but will be calculated exactly once. It's a long calculation so doesn't want to be repeated.
                outputRes = 0.25;
                    #### Calculate binning information
                    #Only do this once because it's computationally expensive but the same for all time steps
                if iCoordMeshes is None: 
                    CCILats = pressure_nc.variables['lat'][:]
                    CCILons = pressure_nc.variables['lon'][:]
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
                timesteps_pressure=len(pressure_time)
                del ilat, ilon, CCILats,CCILons, lat , lon
                
                #### Store data for each day of the month
                pressure_regrid_Vals = np.empty((timesteps_pressure, wind_lat_dimension,wind_lon_dimension), dtype=float);
                # pressure_regrid_ValsErr = np.empty((timesteps_pressure, wind_lat_dimension,wind_lon_dimension), dtype=float);
                # pressure_regrid_ValsCounts = np.empty((timesteps_pressure, wind_lat_dimension,wind_lon_dimension), dtype=float);
                               
                for pressure_timesteps in range(0, timesteps_pressure): 
                    #print(pressure_timesteps)                               
                    pressure_time_slice=pressure_nc.variables['__eo_sp'][pressure_timesteps,:,:]   
                    # no uncertainty data but is needed for function so just use the pressure instead 
                    # and dont use uncertainty.                       
                    pressure_uncertainty_slice=pressure_nc.variables['__eo_sp'][pressure_timesteps,:,:]                         
                    newVals, newCountCount, newValsErr = process_slice(pressure_time_slice,pressure_uncertainty_slice);
                
                    pressure_regrid_Vals[pressure_timesteps,:,:] = newVals;
                    # pressure_regrid_ValsErr[pressure_timesteps,:,:] = newValsErr;
                    # pressure_regrid_ValsCounts[pressure_timesteps,:,:] = newCountCount;
                
                #plt.pcolor(pressure_regrid_Vals[1,:,:])
                             
                pressure_on_wind_grid = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                
                #### get the data of each timestep in wind data
                pressure_time = pressure_nc.variables['time'][:]
                pressure_dates = num2date(pressure_time, pressure_nc.variables['time'].units)
                
                #### loop through the wind timestamps, extract the month and use that to pick which pressure data to use.
                for wind_step in range(0, wind_time_dimension):
                
                    #now find the index of the CLOSEST value and use that
                    abs_deltas_from_target_date = np.absolute(pressure_dates - wind_dates[wind_step])
                    index_of_min_delta_from_target_date = np.argmin(abs_deltas_from_target_date)
                    #closest_date = pressure_dates[index_of_min_delta_from_target_date]

                    pressure_on_wind_grid[wind_step,:,:]=pressure_regrid_Vals[index_of_min_delta_from_target_date,:,:]
                                    
                                                
                #### save pressureitation pre storm output into a netCDF 
                
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
                var = ncout.createVariable("sea_level_pressure", float, ("time","lat", "lon"));
                var.units = "pascals";
                var.long_name = "ERA5 hourly sea level pressure resampled to a 0.25X0.25 degree spatial and hourly temporal resolution";
                var[:] = pressure_on_wind_grid;
                
                ncout.close();   
                
                #### MAXSS ERA5 pressure data pre storm_reference
                
                #pre storm is first 15 days,so 15 days * hourly resolution
                pressure_prestormref=np.nanmean(pressure_on_wind_grid[0:(15*24):1,:,:],axis =(0))
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
                var = ncout.createVariable("sea_level_pressure", float, ("time","lat", "lon"));
                var.units = "pascals";
                var.long_name = "Pre storm (15 days) mean of ERA5 hourly pressure resampled to a 0.25X0.25 degree spatial and hourly temporal resolution";
                var[:] = pressure_on_wind_grid_prestormref;
                
                ncout.close();   
                
                
                #### delete all the variables used during import and saving of pressure      
                del newVals, newValsErr ,pressure_on_wind_grid,ncout,pressure_dates, pressure_lat,pressure_lon,pressure_nc
                del pressure_time,pressure_time_slice,pressure_timesteps,pressure_uncertainty_slice,
                del processedFilePath,timesteps_pressure,abs_deltas_from_target_date,iCoordMeshes,index_of_min_delta_from_target_date
                del  pressure_regrid_Vals, newCountCount
                print("pressure regridded for Storm = "+storm)    
                
                
                ## SCRIPT CHECKED UP TO HERE ##
                
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
                var = ncout.createVariable("land proportion", float, ("latitude", "longitude"));
                var.units = "Percentage";
                var.long_name = "Fraction of grid cell as land (%)";
                var[:] = wind_land_fraction;
                
                ncout.close();   
                
                print("Land fraction for Storm = "+storm)

                

                #### Watson 2020 pCO2 dataset- ADD!!
                
                
                #### Gregor 2021 ETHZ pCO2- ADD!!
                
                
                #### Verification data - pCO2 data monthly SOCAT v4
                downloadedRoot=("verification_data\\SOCATv4");
                downloadedFileTemplate = Template(path.join(downloadedRoot,"2010${MM}01_OCF-CO2-GLO-1M-KRG-CLIM.nc"));
                all_pco2_socatv4 = np.empty((12, 180, 360), dtype=float);
                all_pco2_socatv4[:] = np.nan
                all_conc_co2_air_socatv4 = np.empty((12, 180, 360), dtype=float);
                all_conc_co2_air_socatv4[:] = np.nan
                all_reynolds_socatv4 = np.empty((12, 180, 360), dtype=float);
                all_reynolds_socatv4[:] = np.nan
                
                all_pco2_socatv4_region_subset = np.empty((12, wind_lat_dimension, wind_lon_dimension), dtype=float);
                all_conc_co2_air_socatv4_region_subset = np.empty((12, wind_lat_dimension, wind_lon_dimension), dtype=float);
                all_reynolds_socatv4_region_subset = np.empty((12, wind_lat_dimension, wind_lon_dimension), dtype=float);
                
                    #### resample pco2 to wind spatial grid
                for imonth in range(0, 12):
                    monthStr = format(imonth+1, "02d");
                    downloadFilePath = downloadedFileTemplate.safe_substitute(MM=monthStr);
                    pco2_nc = nc.Dataset(downloadFilePath, 'r')
                    all_pco2_socatv4[:,:] = pco2_nc.variables["fCO2_2010_interpolated_pred"][:,:];
                    all_conc_co2_air_socatv4[:,:] = pco2_nc.variables["vCO2"][:,:];
                    all_reynolds_socatv4[:,:] = pco2_nc.variables["Tcl_2010"][:,:];
                
                    x=all_pco2_socatv4[imonth,:,:]
                    xa=all_conc_co2_air_socatv4[imonth,:,:]
                    xb=all_reynolds_socatv4[imonth,:,:]
                
                    # to do the interpolation, the grid needs to be as column lists, this does that
                    
                    #these lines list lat and long of all points as columns
                    x_coord_range = [i for i in range(0, 360, 1)]
                    y_coord_range = [i for i in range(0, 180, 1)]
                    xy_coord = list(itertools.product(x_coord_range, y_coord_range))
                    
                    #turn 2d data grid of data to column
                    #pco2 sw
                    values = x.flatten(order='F')
                    #conc co2 air
                    values2 = xa.flatten(order='F')
                    # reynolds sst at pco2 seawater
                    values3 = xb.flatten(order='F')
                
                    
                    
                    sample_df = pd.DataFrame()
                
                    
                    #this is first input into interpolation function
                    sample_df['X'] = [xy[0] for xy in xy_coord]
                    sample_df['Y'] = [xy[1] for xy in xy_coord]
                    
                    #this is second input into interpolation function
                    sample_df1 = pd.DataFrame()
                    sample_df1['value'] = values
                    
                    values[values == -999] = 'nan' # or use np.nan
                
                    sample_df2 = pd.DataFrame()
                    sample_df2['value'] = values2
                    values2[values2 == -999] = 'nan' # or use np.nan
                    # get from metadata
                    sample_df3 = pd.DataFrame()
                    sample_df3['value'] = values3
                    values3[values3 == -999] = 'nan' # or use np.nan
                
                    #this is the new grid for the data. global 0.25 degree
                    grid_x, grid_y = np.mgrid[ 0:360:1440j,0:180:720j]
                
                    #grid_z0 = griddata(sample_df, sample_df1, (grid_x,grid_y), method='nearest')
                    grid_z1 = griddata(sample_df, sample_df1, (grid_x,grid_y), method='linear')
                    grid_z1a = griddata(sample_df, sample_df2, (grid_x,grid_y), method='linear')
                    grid_z1b = griddata(sample_df, sample_df3, (grid_x,grid_y), method='linear')
                
                    #grid_z2 = griddata(sample_df, sample_df1, (grid_x,grid_y), method='cubic')
                
                    yyy=grid_z1[:,:,0]
                    yyy=np.rot90(yyy,1)
                    yyy=np.flipud(yyy)
                
                    zzz=grid_z1a[:,:,0]
                    zzz=np.rot90(zzz,1)
                    zzz=np.flipud(zzz)
                
                    nnn=grid_z1b[:,:,0]
                    nnn=np.rot90(nnn,1)
                    nnn=np.flipud(nnn)
                
                    #plt.pcolor(yyy)
                    lon_new_grid = np.arange(-180, 180, 0.25)
                    lat_new_grid = np.arange(-90, 90, 0.25)
                
                    #now work out the subset of the data o extract pco2 for!
                    max_lat_ind = np.where(lat_new_grid == max_lat)[0][0]
                    min_lat_ind = np.where(lat_new_grid == min_lat)[0][0]-1
                
                    max_lon_ind = np.where(lon_new_grid == max_lon)[0][0]
                    min_lon_ind = np.where(lon_new_grid == min_lon)[0][0]-1
                    
                    yyy_subset=yyy[min_lat_ind:max_lat_ind,min_lon_ind:max_lon_ind]
                    zzz_subset=zzz[min_lat_ind:max_lat_ind,min_lon_ind:max_lon_ind]
                    nnn_subset=nnn[min_lat_ind:max_lat_ind,min_lon_ind:max_lon_ind]
                
                    # plt.pcolor(yyy_subset)
                    # plt.pcolor(all_pco2_socatv4_region_subset[1,:,:])
                
                    all_pco2_socatv4_region_subset[imonth,:,:] = yyy_subset;
                    all_conc_co2_air_socatv4_region_subset[imonth,:,:] = zzz_subset;
                    all_reynolds_socatv4_region_subset[imonth,:,:] = nnn_subset;
                
                # make pCO2 matrix the same temporal scale as wind data
                
                #create empty pco2 grid the same size as the wind data.
                pco2_on_wind_grid = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                pco2_on_wind_grid[:] = np.nan
                conc_pco2_air_on_wind_grid = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                pco2_on_wind_grid[:] = np.nan
                reynolds_co2_on_wind_grid = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                pco2_on_wind_grid[:] = np.nan
                
                
                #loop through the wind timestamps, extract the month and use that to pick which pco2 data to use.
                for wind_step in range(0, wind_time_dimension):
                
                    a=wind_dates[wind_step].timetuple()
                    #this is month a[1]
                    #need to index at 0 so -1 month
                    month=a[1]-1
                    
                    pco2_on_wind_grid[wind_step,:,:]=all_pco2_socatv4_region_subset[month,:,:]
                    conc_pco2_air_on_wind_grid[wind_step,:,:]=all_conc_co2_air_socatv4_region_subset[month,:,:]
                    reynolds_co2_on_wind_grid[wind_step,:,:]=all_reynolds_socatv4_region_subset[month,:,:]
                
                
                #### save pco2 output into a netCDF 
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_verification_data_SOCATv4.nc".format(region,year,storm)));
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
                var = ncout.createVariable("pCO2water_mean", float, ("time","lat", "lon"));
                var.units = "uatm";
                var.long_name = "Mean monthly pco2 resampled to hourly on a 0.25X0.25 degree spatial resolution";
                var[:] = pco2_on_wind_grid;
                
                #data variables
                var = ncout.createVariable("xCO2air_mean", float, ("time","lat", "lon"));
                var.units = "ppm";
                var.long_name = "Mean monthly conc co2 in air resampled to hourly on a 0.25X0.25 degree spatial resolution";
                var[:] = conc_pco2_air_on_wind_grid;
                
                #data variables
                var = ncout.createVariable("reynolds_temperature_mean", float, ("time","lat", "lon"));
                var.units = "kelvin";
                var.long_name = "Mean monthly reynolds SST with pco2 seawater resampled to hourly on a 0.25X0.25 degree spatial resolution";
                var[:] = reynolds_co2_on_wind_grid;
                
                ncout.close();   
                
                print("CO2 regridded for Storm = "+storm)
                
                
                
                
                #### Verification data - atm pressure data monthly ECMWF
                downloadedRoot=("verification_data\\air_pressure\\2010");
                downloadedFileTemplate = Template(path.join(downloadedRoot,"2010${MM}_OCF-PRE-GLO-1M-100-ECMWF.nc"));
                all_air_pressure = np.empty((12, 180, 360), dtype=float);
                
                all_air_pressure_region_subset = np.empty((12, wind_lat_dimension, wind_lon_dimension), dtype=float);
                
                
                    #### resample atmospheric pressure to wind spatial grid
                for imonth in range(0, 12):
                    monthStr = format(imonth+1, "02d");
                    downloadFilePath = downloadedFileTemplate.safe_substitute(MM=monthStr);
                    air_pressure_nc = nc.Dataset(downloadFilePath, 'r')
                    all_air_pressure[:,:] = air_pressure_nc.variables["msl_mean"][:,:];
                
                
                    x=all_air_pressure[imonth,:,:]
                
                    # to do the interpolation, the grid needs to be as column lists, this does that
                    
                    #these lines list lat and long of all points as columns
                    x_coord_range = [i for i in range(0, 360, 1)]
                    y_coord_range = [i for i in range(0, 180, 1)]
                    xy_coord = list(itertools.product(x_coord_range, y_coord_range))
                    
                    #turn 2d data grid of data to column
                    #air pressure
                    values = x.flatten(order='F')
                
                    sample_df = pd.DataFrame()
                
                    
                    #this is first input into interpolation function
                    sample_df['X'] = [xy[0] for xy in xy_coord]
                    sample_df['Y'] = [xy[1] for xy in xy_coord]
                    
                    #this is second input into interpolation function
                    sample_df1 = pd.DataFrame()
                    sample_df1['value'] = values
                
                    
                    #this is the new grid for the data. global 0.25 degree
                    grid_x, grid_y = np.mgrid[ 0:360:1440j,0:180:720j]
                
                    #grid_z0 = griddata(sample_df, sample_df1, (grid_x,grid_y), method='nearest')
                    grid_z1 = griddata(sample_df, sample_df1, (grid_x,grid_y), method='linear')
                
                
                    #grid_z2 = griddata(sample_df, sample_df1, (grid_x,grid_y), method='cubic')
                
                    yyy=grid_z1[:,:,0]
                    yyy=np.rot90(yyy,1)
                    yyy=np.flipud(yyy)
                
                
                    #plt.pcolor(yyy)
                    lon_new_grid = np.arange(-180, 180, 0.25)
                    lat_new_grid = np.arange(-90, 90, 0.25)
                
                    #now work out the subset of the data o extract pco2 for!
                    max_lat_ind = np.where(lat_new_grid == max_lat)[0][0]
                    min_lat_ind = np.where(lat_new_grid == min_lat)[0][0]-1
                
                    max_lon_ind = np.where(lon_new_grid == max_lon)[0][0]
                    min_lon_ind = np.where(lon_new_grid == min_lon)[0][0]-1
                    
                    yyy_subset=yyy[min_lat_ind:max_lat_ind,min_lon_ind:max_lon_ind]
                
                    # plt.pcolor(yyy_subset)
                    # plt.pcolor(all_air_pressure_region_subset[1,:,:])
                
                    all_air_pressure_region_subset[imonth,:,:] = yyy_subset;
                
                    #### make pCO2 matrix the same temporal scale as wind data
                
                #create empty pco2 grid the same size as the wind data.
                air_pressure_on_wind_grid = np.empty((wind_time_dimension, wind_lat_dimension, wind_lon_dimension), dtype=float);
                
                
                    #### loop through the wind timestamps, 
                    #extract the month and use that to pick which pressure data to use.
                for wind_step in range(0, wind_time_dimension):
                
                    a=wind_dates[wind_step].timetuple()
                    #this is month a[1]
                    #indexing starts at 0 so minus 1
                    month=a[1]-1
                    
                    air_pressure_on_wind_grid[wind_step,:,:]=all_air_pressure_region_subset[month,:,:]
                
                    #### save air pressure output into a netCDF 
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_verification_data_ECMWF_air_pressure.nc".format(region,year,storm)));
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
                var = ncout.createVariable("sea_level_pressure", float, ("time","lat", "lon"));
                var.units = "pascals";
                var.long_name = "Mean monthly air pressure resampled to hourly on a 0.25X0.25 degree spatial resolution";
                var[:] = air_pressure_on_wind_grid;
                
                ncout.close();   
                print("Air pressure regridded for Storm = "+storm)

                
 
                
                #### world seas mask
                downloadedRoot=("verification_data\\");
                downloadFilePath = path.join(downloadedRoot,"World_Seas-final-complete_IGA.nc");
                world_seas = np.empty((180, 360), dtype=float);
                world_seas_nc = nc.Dataset(downloadFilePath, 'r')
                world_seas[:,:] = world_seas_nc.variables["sea-mask"][:,:];
                x=world_seas[:,:]

                
                #these lines list lat and long of all points as columns
                x_coord_range = [i for i in range(0, 360, 1)]
                y_coord_range = [i for i in range(0, 180, 1)]
                xy_coord = list(itertools.product(x_coord_range, y_coord_range))
                
                #turn 2d data grid of data to column
                values = x.flatten(order='F')
                 
                sample_df = pd.DataFrame()

                #this is first input into interpolation function
                sample_df['X'] = [xy[0] for xy in xy_coord]
                sample_df['Y'] = [xy[1] for xy in xy_coord]
                
                #this is second input into interpolation function
                sample_df1 = pd.DataFrame()
                sample_df1['value'] = values
                
                values[values == -999] = 'nan' # or use np.nan
        
                #this is the new grid for the data. global 0.25 degree
                grid_x, grid_y = np.mgrid[ 0:360:1440j,0:180:720j]
            
                #grid_z0 = griddata(sample_df, sample_df1, (grid_x,grid_y), method='nearest')
                #grid_z1 = griddata(sample_df, sample_df1, (grid_x,grid_y), method='linear')
                grid_z1 = griddata(sample_df, sample_df1, (grid_x,grid_y), method='nearest')

                #grid_z2 = griddata(sample_df, sample_df1, (grid_x,grid_y), method='cubic')
            
                yyy=grid_z1[:,:,0]
                yyy=np.rot90(yyy,1)
                #yyy=np.flipud(yyy)
            
                #plt.pcolor(yyy)
                lon_new_grid = np.arange(-180, 180, 0.25)
                lat_new_grid = np.arange(-90, 90, 0.25)
            
                #now work out the subset of the data to extract for!
                max_lat_ind = np.where(lat_new_grid == max_lat)[0][0]
                min_lat_ind = np.where(lat_new_grid == min_lat)[0][0]-1
            
                max_lon_ind = np.where(lon_new_grid == max_lon)[0][0]
                min_lon_ind = np.where(lon_new_grid == min_lon)[0][0]-1
                
                yyy_subset=yyy[min_lat_ind:max_lat_ind,min_lon_ind:max_lon_ind]

            

                #### save mask output into a netCDF 
                processedFilePath = (path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\Resampled_for_fluxengine_world_ocean_mask.nc".format(region,year,storm)));
                ncout = Dataset(processedFilePath, 'w');
                # create dataset and provide dimensions
                
                ncout.createDimension("lat", wind_lat_dimension);
                ncout.createDimension("lon", wind_lon_dimension);
                
                #dimension variables
                var = ncout.createVariable("lat", float, ("lat",));
                var.units = "lat (degrees North)";
                var[:] = wind_lat;
                
                var = ncout.createVariable("lon", float, ("lon",));
                var.units = "lon (degrees East)";
                var[:] = wind_lon;
                

                #data variables
                var = ncout.createVariable("seamask", float, ("lat", "lon"));
                var.units = "unitless";
                var.long_name = "World seas mask nearest interpolation on to  a 0.25X0.25 degree spatial resolution";
                var[:] = yyy_subset;
                
                ncout.close();   
                
                print("mask regridded for Storm = "+storm)

                
                
                

# -*- coding: utf-8 -*-
"""
Script for looping through the MAXSS storm-Atlas and running fluxengine

@author: Richard Sims
"""

    #### all functions and scripts imported here


#fluxengine is installed here. 
#C:\Users\rps207\Documents\Python\2021-Anaconda_install\Lib\site-packages\fluxengine


#pacakges
import os
from os import path;
from fluxengine.core.fe_setup_tools import run_fluxengine, get_fluxengine_root;
from fluxengine.tools.lib_ofluxghg_flux_budgets import run_flux_budgets;
from fluxengine.tools.lib_compare_net_budgets import read_global_core_budgets, calc_net_budget_percentages;
from glob import glob
from pathlib import Path
import shutil
import netCDF4 as nc
from netCDF4 import date2num, num2date, Dataset
from argparse import Namespace;
import numpy as np
from datetime import datetime, timedelta;
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from pyproj import Geod # use pyproj as it is documented code

import inspect;
Month_Fmt = mdates.DateFormatter('%b %d')


def get_datetimes(secondsSince1970):
    base = datetime(1970,1,1);
    return np.array([base+timedelta(seconds=int(t)) for t in secondsSince1970]);

def make_configuration_file(storm_dir_relative,timestepsinfile,region,year,storm,run_name): # 

    # if it doesn't exist make directory for config files 
    config_folder_Path=path.join("output\\configs\\{2}\\maxss\\storm-atlas\\ibtracs\\{0}\\{1}".format(region,year,run_name))
    if not os.path.exists(config_folder_Path):
        os.makedirs(config_folder_Path)

    #copy configuration file template
    configfiletemplate="E:/MAXSS_Wilson/MAXSS_configuration_file_template.conf"
    
    configfilenew = os.path.join(config_folder_Path, f"{storm}.conf")
    shutil.copy(configfiletemplate,configfilenew)


    #there is now a configuration file specific for this storm
    #time to edit that config file with the correct paths!

    # Read in the file
    with open(configfilenew, 'r') as file :
      filedata = file.read()
    
    # Replace the target strings in configuration file template with new paths
    # or info that changes between storms
    

    # Wind paths
    if run_name=="MAXSS_RUN" or run_name=="WIND_RUN":
        filedata = filedata.replace('windu10_path = windu10path.nc', 'windu10_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_MAXSS_L4_windspeed.nc')
        filedata = filedata.replace('windu10_temporalChunking = numberoftimesteps','windu10_temporalChunking ='+str(timestepsinfile))
    else:
        filedata = filedata.replace('windu10_path = windu10path.nc', 'windu10_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_MAXSS_L4_windspeed_pre_storm_reference.nc')
        filedata = filedata.replace('windu10_temporalChunking = numberoftimesteps','windu10_temporalChunking ='+str(timestepsinfile))
    

    # Wind moment 2 path
    if run_name=="MAXSS_RUN" or run_name=="WIND_RUN":
        filedata = filedata.replace('windu10_moment2_path = windmoment2path.nc', 'windu10_moment2_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_MAXSS_L4_windspeed.nc')
        filedata = filedata.replace('windu10_moment2_temporalChunking = numberoftimesteps','windu10_moment2_temporalChunking ='+str(timestepsinfile))
    else:
        filedata = filedata.replace('windu10_moment2_path = windmoment2path.nc', 'windu10_moment2_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_MAXSS_L4_windspeed_pre_storm_reference.nc')
        filedata = filedata.replace('windu10_moment2_temporalChunking = numberoftimesteps','windu10_moment2_temporalChunking ='+str(timestepsinfile))


    # Wind moment 3 path
    if run_name=="MAXSS_RUN" or run_name=="WIND_RUN":
        filedata = filedata.replace('windu10_moment3_path = windmoment3path.nc', 'windu10_moment3_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_MAXSS_L4_windspeed.nc')
        filedata = filedata.replace('windu10_moment3_temporalChunking = numberoftimesteps','windu10_moment3_temporalChunking ='+str(timestepsinfile))
    else:
        filedata = filedata.replace('windu10_moment3_path = windmoment3path.nc', 'windu10_moment3_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_MAXSS_L4_windspeed_pre_storm_reference.nc')
        filedata = filedata.replace('windu10_moment3_temporalChunking = numberoftimesteps','windu10_moment3_temporalChunking ='+str(timestepsinfile))
    
    # Wind moment 3.7 path
    if run_name=="MAXSS_RUN" or run_name=="WIND_RUN":
        filedata = filedata.replace('windu10_momentthreeseven_path = windmomentthreesevenpath.nc', 'windu10_momentthreeseven_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_MAXSS_L4_windspeed.nc')
        filedata = filedata.replace('windu10_momentthreeseven_temporalChunking = numberoftimesteps','windu10_momentthreeseven_temporalChunking ='+str(timestepsinfile))
    else:
        filedata = filedata.replace('windu10_momentthreeseven_path = windmomentthreesevenpath.nc', 'windu10_momentthreeseven_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_MAXSS_L4_windspeed_pre_storm_reference.nc')
        filedata = filedata.replace('windu10_momentthreeseven_temporalChunking = numberoftimesteps','windu10_momentthreeseven_temporalChunking ='+str(timestepsinfile))
        

    # SST paths
    if run_name=="MAXSS_RUN" or run_name=="SST_RUN":
        filedata = filedata.replace('sstfnd_path = fndpath.nc', 'sstfnd_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_MAXSS_ESACCI_SST.nc')
        filedata = filedata.replace('sstfnd_temporalChunking = numberoftimesteps','sstfnd_temporalChunking ='+str(timestepsinfile))
    else:
        filedata = filedata.replace('sstfnd_path = fndpath.nc', 'sstfnd_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_MAXSS_ESACCI_SST_pre_storm_reference.nc')
        filedata = filedata.replace('sstfnd_temporalChunking = numberoftimesteps','sstfnd_temporalChunking ='+str(timestepsinfile))
    
    
    # salinity path
    if run_name=="MAXSS_RUN" or run_name=="SSS_RUN":
        filedata = filedata.replace('salinity_path = salinitypath.nc', 'salinity_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_MAXSS_ESACCI_SSS.nc')
        filedata = filedata.replace('salinity_temporalChunking = numberoftimesteps','salinity_temporalChunking ='+str(timestepsinfile))
    else:
        filedata = filedata.replace('salinity_path = salinitypath.nc', 'salinity_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_MAXSS_ESACCI_SSS_pre_storm_reference.nc')
        filedata = filedata.replace('salinity_temporalChunking = numberoftimesteps','salinity_temporalChunking ='+str(timestepsinfile))

    # Pressure path
    if run_name=="MAXSS_RUN" or run_name=="PRESSURE_RUN":
        filedata = filedata.replace('pressure_path = pressurepath.nc', 'pressure_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_MAXSS_ERA5_pressure.nc')
        filedata = filedata.replace('pressure_temporalChunking = numberoftimesteps','pressure_temporalChunking ='+str(timestepsinfile))
        #matching the variable name to that shown in Panoply
        filedata = filedata.replace('pressure_prod = Sea level pressure', 'pressure_prod = sea_level_pressure')
    else:
        filedata = filedata.replace('pressure_path = pressurepath.nc', 'pressure_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_MAXSS_ERA5_pressure_pre_storm_reference.nc')
        filedata = filedata.replace('pressure_temporalChunking = numberoftimesteps','pressure_temporalChunking ='+str(timestepsinfile))
        #matching the variable name to that shown in Panoply
        filedata = filedata.replace('pressure_prod = Sea level pressure', 'pressure_prod = sea_level_pressure')
    # Precipitation path
    if run_name=="MAXSS_RUN" or run_name=="PRECIPITATION_RUN":
        filedata = filedata.replace('rain_path = precipitationpath.nc', 'rain_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_MAXSS_ERA5_precipitation.nc')
        filedata = filedata.replace('rain_temporalChunking = numberoftimesteps','rain_temporalChunking ='+str(timestepsinfile))
    else:
        filedata = filedata.replace('rain_path = precipitationpath.nc', 'rain_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_MAXSS_ERA5_precipitation_pre_storm_reference.nc')
        filedata = filedata.replace('rain_temporalChunking = numberoftimesteps','rain_temporalChunking ='+str(timestepsinfile))


    # xCO2 air path
    filedata = filedata.replace('vgas_air_path = co2airmixingrationpath.nc', 'vgas_air_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_Ford_et_al_pco2.nc')
    filedata = filedata.replace('vgas_air_temporalChunking = numberoftimesteps','vgas_air_temporalChunking ='+str(timestepsinfile))
    #update variable name
    filedata = filedata.replace('vgas_air_prod = xCO2air_mean', 'vgas_air_prod = V_gas')
    
    
    # pco2 path
    filedata = filedata.replace('pgas_sw_path = pco2seawaterpath.nc', 'pgas_sw_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_Ford_et_al_pco2.nc')
    filedata = filedata.replace('pgas_sw_temporalChunking = numberoftimesteps','pgas_sw_temporalChunking ='+str(timestepsinfile))

    # pco2 SST path
    filedata = filedata.replace('pco2_sst_path = pco2swpath.nc', 'pco2_sst_path ='+storm_dir_relative+r'\Resampled_for_fluxengine_Ford_et_al_pco2.nc')
    filedata = filedata.replace('pco2_sst_temporalChunking = numberoftimesteps','pco2_sst_temporalChunking ='+str(timestepsinfile))


    #output directory
    filedata = filedata.replace('output_dir = output/', 'output_dir ='+"output\\{3}\\maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}".format(region,year,storm,run_name))
    #output file format
    filedata = filedata.replace('output_file = MAXSS_DOY_<DDD>_DATE_<YYYY>_<MM>_<DD>.nc', 'output_file = MAXSS_'+storm+'_DOY_<DDD>_DATE_<YYYY>_<MM>_<DD>.nc')

    # Write the file out again
    with open(configfilenew, 'w') as file:
      file.write(filedata)
      
      
    return configfilenew;

def get_spatially_integrated_flux(fe,region,year,storm,run_name,wind_time, wind_storm_land_fraction):
    
    #### Calculate grid cell areas
    rez=0.25 #spatial resolution of grid
    storm_lon=fe.longitude_data # latitude of fluxengine output
    storm_lat=fe.latitude_data # latitude of fluxengine output
    areagrid=0*fe.latitude_grid # lat and longitude grid size, set to 0.

    geod = Geod(ellps='WGS84')        # use PYPROJ and WGS84 - both well documented                      

    lat_counter=0
    for lat in storm_lat: #loop through latitude
        lon_counter=0
        for lon in storm_lon: #loop through longitude
            #define the latitude and longitude coordinates of the four points 
            #for which the lat and lon is the centre. A fifth value is needed
            #to complete the shape.
            lats = [lat-rez/2, lat-rez/2, lat+rez/2, lat+rez/2, lat-rez/2]
            lons = [lon-rez/2, lon+rez/2, lon+rez/2, lon-rez/2,lon-rez/2]
               
            poly_area, poly_perimeter = geod.polygon_area_perimeter(lons, lats)
            #print("area: {} , perimeter: {}".format(poly_area, poly_perimeter))
            areagrid[lat_counter,lon_counter]=poly_area # add the lat and lon values to matrix
            lon_counter=lon_counter+1
        lat_counter=lat_counter+1
            
    #### Calculate Integrated flux
    print("\n\nNow calculating flux budgets for Region={0} year={1} Storm={2}".format(region,year,storm));
    #Flux is m^2 per day
    #Areas are all in m^2
    #Timestep needs to be scaled by data resolution 
    #Equation is then Flux*Area*(timestep in hours/24)
    #The total flux is then the sum of these fluxes
    
    # Get all files and directories ending with .nc
    Fe_oututfile_list = glob(os.path.join(fe.runParams.output_dir, "*.nc"))
    
    #get the temporal resolution from the fluxengine
    #turn it into hours
    timestep=fe.runParams.temporal_resolution
    seconds = timestep.total_seconds()
    fe_temporal_hours = seconds // 3600

    #create variables the size of the hourly timesteps
    Storm_flux_hourly=[]
    Storm_time_hourly=[]

    #loop through netCDF files
    for fluxfile_number in range(0, len(Fe_oututfile_list)):
        #Load using with so that .nc file automatically closed after use
        with nc.Dataset(Fe_oututfile_list[fluxfile_number]) as flux_nc:
                     
            #load in the flux data
            Flux_data=flux_nc.variables['OF'][:,:] #DOES THIS NEED CHANGING? [:] (Appeared not to impact results)
            #now scale the fluxes
            Scaledfluxes=Flux_data*areagrid*(1-wind_storm_land_fraction)*(fe_temporal_hours/24)# g C hr-1 per unit area of the grid cell
            #now sum the fluxes over the whole region (Dan changed to nansum to allow calculation where nan values present)
            Integrated_regional_flux=np.nansum(Scaledfluxes,axis =(1,2))#sum over spatial dimension
            Time_data=flux_nc.variables['time'][:]
            
            #appendthe integrated flux and time to a combined Matrix.
            Storm_flux_hourly.append(Integrated_regional_flux)
            Storm_time_hourly.append(Time_data)
    
    #these are lists- want them as 1d arrays
    #first convert to 2d array
    arr = np.array(Storm_flux_hourly)
    time_arr = np.array(Storm_time_hourly)
    #get size of matrix 
    arr_size=(np.shape(arr))
    #second get dimensions and make 1d array
    arr2=arr.reshape(arr_size[0]*arr_size[1],1)
    #convert from grams to Tg
    unit_factor=(1000000000000)
    #Integrated flux in Tg every hour
    Hourly_regional_flux_Tg=arr2/unit_factor
    #reshape time array
    time_arr2=time_arr.reshape(arr_size[0]*arr_size[1],1)
    time_for_plotting=get_datetimes(time_arr2)

    #total flux in region across all timesteps
    Regiona_flux_Tg=np.nansum(Hourly_regional_flux_Tg) #Dan changed this to nansum to ensure data with nans was summed
    
    
    #### Save these to a folder
    # if it doesn't exist make directory for Integrated flux files 
    Int_flux_folder_Path=path.join("output\\Spatially_integrated_fluxes\\maxss\\storm-atlas\\ibtracs\\{0}\\{1}".format(region,year,run_name))
    if not os.path.exists(Int_flux_folder_Path):
        os.makedirs(Int_flux_folder_Path)

    processedFilePath = os.path.join(Int_flux_folder_Path, f"{storm}_{run_name}.nc")

    ncout = Dataset(processedFilePath, 'w');
    
        #### provide dimensions
    ncout.createDimension("time", len(time_arr2));

    var = ncout.createVariable("time", int, ("time",));
    var.long_name = "Time";
    var.units = "seconds since 1970-01-01 00:00:00";
    var[:] = time_arr2
    
    #data variables
    var = ncout.createVariable("Hourly_flux", float, ("time",));
    var.units = "Tg C hr-1";
    var.long_name = "Total hourly flux in region";
    var[:] = Hourly_regional_flux_Tg;
    
    #data variables
    var = ncout.createVariable("Total_flux", float);
    var.units = "Tg C";
    var.long_name = "Total flux in region for period of the storm";
    var[:] = Regiona_flux_Tg;
    
    ncout.close();  

    return Hourly_regional_flux_Tg,time_for_plotting;


if __name__ == "__main__":
    verbose=True
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
                storm_dir_relative = path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}".format(region,year,storm));
        
                # need to get the identifier from the storm name as it is used in .nc file string
                storm_id=storm.split('_')[1]
                
                if region=="north-atlantic":
                    region_id="NA"
                elif region=="east-pacific":
                    region_id="EP"
                elif region=="north-indian":
                    region_id="NI"
                elif region=="south-atlantic":
                    region_id="SA"
                elif region=="south-indian":
                    region_id="SI"
                elif region=="south-pacific":
                    region_id="SP"
                elif region=="west-pacific":
                    region_id="WP"  
                    
                    
                #### Need to add temporal chunking to config file
                #- so need to open wind data to get that
                # Use 'with' to ensure the wind file closes immediately after extracting metadata
                with nc.Dataset(path.join("maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}\\MAXSS_HIST_TC_{3}_{1}_{4}_MAXSS_HIST_TC_L4.nc".format(region,year,storm,region_id,storm_id))) as winds_nc:
                    wind_northward = winds_nc.variables['__eo_northward_wind'][:]
                    # need land fraction mask from wind data
                    wind_storm_land_fraction = winds_nc.variables['__eo_land_fraction'][0]
                    wind_time_dimension=len(wind_northward)
                    timestepsinfile=wind_time_dimension
                    
                    #### get start and end times - to be added to run call
                    wind_time = winds_nc.variables['time'][:]
                    wind_dates = num2date(wind_time, winds_nc.variables['time'].units)

                run_startime=wind_dates[0].strftime("%Y-%m-%d %H:%M")#
                run_endtime=wind_dates[-1].strftime("%Y-%m-%d %H:%M")#
          
                #run_endtime=wind_dates[24].strftime("%Y-%m-%d %H:%M")# # TO ONLY RUN FOR one day
                

                #### Run flux engine for 'MAXSS run'
                run_name="MAXSS_RUN"
                # create custom config file for this storm
                # call custom function which copies file template and makes edits
                configFilePath_MAXSS_RUN=make_configuration_file(storm_dir_relative,timestepsinfile,region,year,storm,run_name)
                print("Running FluxEngine for Region={0} year={1} Storm={2}".format(region,year,storm));
                runStatus, fe_MAXSS_RUN = run_fluxengine(configFilePath_MAXSS_RUN,run_startime,run_endtime,processLayersOff=True, verbose=False);
                #call function to get sum of hourly fluxes scaled by area.
                Hourlyflux_MAXSS_RUN,Hourlyfluxdate_MAXSS_RUN=get_spatially_integrated_flux(fe_MAXSS_RUN,region,year,storm,run_name,wind_time, wind_storm_land_fraction)

                #### Run flux engine for "REF run"
                run_name="REF_RUN"
                # create custom config file for this storm
                # call custom function which copies file template and makes edits
                configFilePath_REF_RUN=make_configuration_file(storm_dir_relative,timestepsinfile,region,year,storm,run_name)
                print("Running FluxEngine for Region={0} year={1} Storm={2}".format(region,year,storm));
                runStatus, fe_REF_RUN = run_fluxengine(configFilePath_REF_RUN,run_startime,run_endtime,processLayersOff=True, verbose=False);
                #call function to get sum of hourly fluxes scaled by area.
                Hourlyflux_REF_RUN,Hourlyfluxdate_REF_RUN=get_spatially_integrated_flux(fe_REF_RUN,region,year,storm,run_name,wind_time, wind_storm_land_fraction)

                
                #### Run flux engine for "WIND run"
                run_name="WIND_RUN"
                # create custom config file for this storm
                # call custom function which copies file template and makes edits
                configFilePath_WIND_RUN=make_configuration_file(storm_dir_relative,timestepsinfile,region,year,storm,run_name)
                print("Running FluxEngine for Region={0} year={1} Storm={2}".format(region,year,storm));
                runStatus, fe_WIND_RUN = run_fluxengine(configFilePath_WIND_RUN,run_startime,run_endtime,processLayersOff=True, verbose=False);
                #call function to get sum of hourly fluxes scaled by area.
                Hourlyflux_WIND_RUN,Hourlyfluxdate_WIND_RUN=get_spatially_integrated_flux(fe_WIND_RUN,region,year,storm,run_name,wind_time, wind_storm_land_fraction)

                #### Run flux engine for "SST run"
                run_name="SST_RUN"
                # create custom config file for this storm
                # call custom function which copies file template and makes edits
                configFilePath_SST_RUN=make_configuration_file(storm_dir_relative,timestepsinfile,region,year,storm,run_name)
                print("Running FluxEngine for Region={0} year={1} Storm={2}".format(region,year,storm));
                runStatus, fe_SST_RUN = run_fluxengine(configFilePath_SST_RUN,run_startime,run_endtime,processLayersOff=True, verbose=False);
                #call function to get sum of hourly fluxes scaled by area.
                Hourlyflux_SST_RUN,Hourlyfluxdate_SST_RUN=get_spatially_integrated_flux(fe_SST_RUN,region,year,storm,run_name,wind_time)

                #### Run flux engine for "SSS run"
                run_name="SSS_RUN"
                # create custom config file for this storm
                # call custom function which copies file template and makes edits
                configFilePath_SSS_RUN=make_configuration_file(storm_dir_relative,timestepsinfile,region,year,storm,run_name)
                print("Running FluxEngine for Region={0} year={1} Storm={2}".format(region,year,storm));
                runStatus, fe_SSS_RUN = run_fluxengine(configFilePath_SSS_RUN,run_startime,run_endtime,processLayersOff=True, verbose=False);
                #call function to get sum of hourly fluxes scaled by area.
                Hourlyflux_SSS_RUN,Hourlyfluxdate_SSS_RUN=get_spatially_integrated_flux(fe_SSS_RUN,region,year,storm,run_name,wind_time, wind_storm_land_fraction)

            
                #### Run flux engine for "PRESSURE run"
                run_name="PRESSURE_RUN"
                # create custom config file for this storm
                # call custom function which copies file template and makes edits
                configFilePath_PRESSURE_RUN=make_configuration_file(storm_dir_relative,timestepsinfile,region,year,storm,run_name)
                print("Running FluxEngine for Region={0} year={1} Storm={2}".format(region,year,storm));
                runStatus, fe_PRESSURE_RUN = run_fluxengine(configFilePath_PRESSURE_RUN,run_startime,run_endtime,processLayersOff=True, verbose=False);
                #call function to get sum of hourly fluxes scaled by area.
                Hourlyflux_PRESSURE_RUN,Hourlyfluxdate_PRESSURE_RUN=get_spatially_integrated_flux(fe_PRESSURE_RUN,region,year,storm,run_name,wind_time, wind_storm_land_fraction)

                
                #### Run flux engine for "PRECIPITATION run"
                run_name="PRECIPITATION_RUN"
                # create custom config file for this storm
                # call custom function which copies file template and makes edits
                configFilePath_PRECIPITATION_RUN=make_configuration_file(storm_dir_relative,timestepsinfile,region,year,storm,run_name)
                print("Running FluxEngine for Region={0} year={1} Storm={2}".format(region,year,storm));
                runStatus, fe_PRECIPITATION_RUN = run_fluxengine(configFilePath_PRECIPITATION_RUN,run_startime,run_endtime,processLayersOff=True, verbose=False);
                #call function to get sum of hourly fluxes scaled by area.
                Hourlyflux_PRECIPITATION_RUN,Hourlyfluxdate_PRECIPITATION_RUNN=get_spatially_integrated_flux(fe_PRECIPITATION_RUN,region,year,storm,run_name,wind_time, wind_storm_land_fraction)


                # Add to storm counter when everything is done.
                storm_counter=storm_counter+1
                



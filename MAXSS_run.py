# -*- coding: utf-8 -*-
"""
Script for looping through the MAXSS storm-Atlas and running fluxengine

@author: Richard Sims
"""
#Main script for running the FluxEngine TC simulations

#fluxengine is installed here.
#C:\Users\rps207\Documents\Python\2021-Anaconda_install\Lib\site-packages\fluxengine

# test git working correctly #

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
import argparse #so we can run this code from bash script

import inspect;
Month_Fmt = mdates.DateFormatter('%b %d')

variables_to_exclude = [
        "pressure", "sstfnd", "sstfndC", "sstfnd_count", "sstfnd_stddev",
        "pco2_sst", "rain", "sstskin", "sstskin_count", "sstskin_stddev",
        "vgas_air", "windu10_moment2", "windu10_moment3",
        "windu10_momentthreeseven", "windu10"]


def get_datetimes(secondsSince1970):
    base = datetime(1970,1,1);
    return np.array([base+timedelta(seconds=int(t)) for t in secondsSince1970]);

def make_configuration_file(storm_dir_relative,timestepsinfile,region,year,storm,run_name,time_chunk_val,configfiletemplate="E:/MAXSS_Wilson/MAXSS_configuration_file_template.conf"): #
    #DJF 09/05/2026: Added the configfiletemplate as a function input so it can be modified by the main function

    # if it doesn't exist make directory for config files
    config_folder_path = path.join("output", "configs", run_name, "maxss","storm-atlas","tropical","ibtracs", region, year)
    
    if not os.path.exists(config_folder_path):
        os.makedirs(config_folder_path)

    #copy configuration file template
    # configfiletemplate="E:/MAXSS_Wilson/MAXSS_configuration_file_template.conf" # DJF 09/05/2026: Added as function input

    configfilenew = os.path.join(config_folder_path, f"{storm}.conf")
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
        filedata = filedata.replace('windu10_path = windu10path.nc', 'windu10_path ='+ storm_dir_relative + '/Resampled_for_fluxengine_MAXSS_L4_windspeed.nc')
        filedata = filedata.replace('windu10_temporalChunking = numberoftimesteps','windu10_temporalChunking ='+str(timestepsinfile))
    else:
        filedata = filedata.replace('windu10_path = windu10path.nc', 'windu10_path ='+ storm_dir_relative + '/Resampled_for_fluxengine_MAXSS_L4_windspeed_pre_storm_reference.nc')
        filedata = filedata.replace('windu10_temporalChunking = numberoftimesteps','windu10_temporalChunking ='+str(timestepsinfile))

    # Wind moment 2 path
    if run_name=="MAXSS_RUN" or run_name=="WIND_RUN":
        filedata = filedata.replace('windu10_moment2_path = windmoment2path.nc', 'windu10_moment2_path ='+storm_dir_relative+ '/Resampled_for_fluxengine_MAXSS_L4_windspeed.nc')
        filedata = filedata.replace('windu10_moment2_temporalChunking = numberoftimesteps','windu10_moment2_temporalChunking ='+str(timestepsinfile))
    else:
        filedata = filedata.replace('windu10_moment2_path = windmoment2path.nc', 'windu10_moment2_path ='+storm_dir_relative+ '/Resampled_for_fluxengine_MAXSS_L4_windspeed_pre_storm_reference.nc')
        filedata = filedata.replace('windu10_moment2_temporalChunking = numberoftimesteps','windu10_moment2_temporalChunking ='+str(timestepsinfile))

    # SST (with and without gradients) paths
    if run_name=="MAXSS_RUN" or run_name=="SST_WITH_GRADIENTS_RUN" or run_name=="SST_NO_GRADIENTS_RUN":
        # SST variable set up
        filedata = filedata.replace('sstfnd_path = fndpath.nc', 'sstfnd_path ='+storm_dir_relative+ '/Resampled_for_fluxengine_MAXSS_ESACCI_SST.nc')
        filedata = filedata.replace('sstfnd_temporalChunking = numberoftimesteps','sstfnd_temporalChunking ='+str(timestepsinfile))

        # pco2 variable set up
        filedata = filedata.replace('pgas_sw_path = pco2seawaterpath.nc', 'pgas_sw_path ='+storm_dir_relative+ '/Resampled_for_fluxengine_Ford_et_al_pco2.nc')
        filedata = filedata.replace('pgas_sw_temporalChunking = numberoftimesteps','pgas_sw_temporalChunking ='+str(timestepsinfile))

        # pco2 SST variable set up
        filedata = filedata.replace('pco2_sst_path = pco2swpath.nc', 'pco2_sst_path ='+storm_dir_relative+ '/Resampled_for_fluxengine_Ford_et_al_pco2.nc')
        filedata = filedata.replace('pco2_sst_temporalChunking = numberoftimesteps','pco2_sst_temporalChunking ='+str(timestepsinfile))

    else:
        # SST variable set up
        filedata = filedata.replace('sstfnd_path = fndpath.nc', 'sstfnd_path ='+storm_dir_relative+ '/Resampled_for_fluxengine_MAXSS_ESACCI_SST_pre_storm_reference.nc')
        filedata = filedata.replace('sstfnd_temporalChunking = numberoftimesteps','sstfnd_temporalChunking ='+str(timestepsinfile))

        # pco2 variable set up
        filedata = filedata.replace('pgas_sw_path = pco2seawaterpath.nc', 'pgas_sw_path ='+storm_dir_relative+ '/Resampled_for_fluxengine_Ford_et_al_pco2_pre_storm_reference.nc')
        filedata = filedata.replace('pgas_sw_temporalChunking = numberoftimesteps','pgas_sw_temporalChunking ='+str(timestepsinfile))

        # pco2 SST variable set up
        filedata = filedata.replace('pco2_sst_path = pco2swpath.nc', 'pco2_sst_path ='+storm_dir_relative+ '/Resampled_for_fluxengine_Ford_et_al_pco2_pre_storm_reference.nc')
        filedata = filedata.replace('pco2_sst_temporalChunking = numberoftimesteps','pco2_sst_temporalChunking ='+str(timestepsinfile))

    # SST no gradients (turn off gradients)
    #(Default in the config file is that gradients turned on, so we need to switch off for this run)
    if run_name=="SST_NO_GRADIENTS_RUN" :
        filedata = filedata.replace('sst_gradients = yes', 'sst_gradients = no')
        filedata = filedata.replace('cool_skin_difference = 0.17', 'cool_skin_difference = 0')
        filedata = filedata.replace('saline_skin_value = 0.1', 'saline_skin_value = 0')

    if run_name=="MAXSS_RUN" or run_name=="V_GAS_RUN":
        # V_gas path
        filedata = filedata.replace('vgas_air_path = co2airmixingrationpath.nc', 'vgas_air_path ='+storm_dir_relative+ '/Resampled_for_fluxengine_Ford_et_al_pco2.nc')
        filedata = filedata.replace('vgas_air_temporalChunking = numberoftimesteps','vgas_air_temporalChunking ='+str(timestepsinfile))
        #update variable name
        filedata = filedata.replace('vgas_air_prod = xCO2air_mean', 'vgas_air_prod = V_gas')
    else:
        filedata = filedata.replace('vgas_air_path = co2airmixingrationpath.nc', 'vgas_air_path ='+storm_dir_relative+ '/Resampled_for_fluxengine_Ford_et_al_pco2_pre_storm_reference.nc')
        filedata = filedata.replace('vgas_air_temporalChunking = numberoftimesteps','vgas_air_temporalChunking ='+str(timestepsinfile))
        #update variable name
        filedata = filedata.replace('vgas_air_prod = xCO2air_mean', 'vgas_air_prod = V_gas')

    # salinity path
    if run_name=="MAXSS_RUN" or run_name=="SSS_RUN":
        filedata = filedata.replace('salinity_path = salinitypath.nc', 'salinity_path ='+storm_dir_relative+ '/Resampled_for_fluxengine_MAXSS_ESACCI_SSS.nc')
        filedata = filedata.replace('salinity_temporalChunking = numberoftimesteps','salinity_temporalChunking ='+str(timestepsinfile))
    else:
        filedata = filedata.replace('salinity_path = salinitypath.nc', 'salinity_path ='+storm_dir_relative+ '/Resampled_for_fluxengine_MAXSS_ESACCI_SSS_pre_storm_reference.nc')
        filedata = filedata.replace('salinity_temporalChunking = numberoftimesteps','salinity_temporalChunking ='+str(timestepsinfile))

    # Pressure path
    if run_name=="MAXSS_RUN" or run_name=="PRESSURE_RUN":
        filedata = filedata.replace('pressure_path = pressurepath.nc', 'pressure_path ='+storm_dir_relative+ '/Resampled_for_fluxengine_MAXSS_ERA5_pressure.nc')
        filedata = filedata.replace('pressure_temporalChunking = numberoftimesteps','pressure_temporalChunking ='+str(timestepsinfile))
        #matching the variable name to that shown in Panoply
        filedata = filedata.replace('pressure_prod = Sea level pressure', 'pressure_prod = sea_level_pressure')
    else:
        filedata = filedata.replace('pressure_path = pressurepath.nc', 'pressure_path ='+storm_dir_relative+ '/Resampled_for_fluxengine_MAXSS_ERA5_pressure_pre_storm_reference.nc')
        filedata = filedata.replace('pressure_temporalChunking = numberoftimesteps','pressure_temporalChunking ='+str(timestepsinfile))
        #matching the variable name to that shown in Panoply
        filedata = filedata.replace('pressure_prod = Sea level pressure', 'pressure_prod = sea_level_pressure')
    # Precipitation path
    if run_name=="MAXSS_RUN" or run_name=="PRECIPITATION_RUN":
        filedata = filedata.replace('rain_path = precipitationpath.nc', 'rain_path ='+storm_dir_relative+ '/Resampled_for_fluxengine_MAXSS_ERA5_precipitation.nc')
        filedata = filedata.replace('rain_temporalChunking = numberoftimesteps','rain_temporalChunking ='+str(timestepsinfile))
    else:
        filedata = filedata.replace('rain_path = precipitationpath.nc', 'rain_path ='+storm_dir_relative+ '/Resampled_for_fluxengine_MAXSS_ERA5_precipitation_pre_storm_reference.nc')
        filedata = filedata.replace('rain_temporalChunking = numberoftimesteps','rain_temporalChunking ='+str(timestepsinfile))

    #output directory
    filedata = filedata.replace('output_dir = output/', 'output_dir = ' + "output/{3}/maxss/storm-atlas/tropical/ibtracs/{0}/{1}/{2}".format(region, year, storm, run_name))
    
    #output file format
    filedata = filedata.replace('output_file = MAXSS_DOY_<DDD>_DATE_<YYYY>_<MM>_<DD>.nc', 'output_file = MAXSS_'+storm+'_DOY_<DDD>_DATE_<YYYY>_<MM>_<DD>.nc')

    # Exclude variables that are not needed from the netcdf outputs
    filedata += "\nexclude_outputs = " + ",".join(variables_to_exclude) + "\n"

    #Set the mask used to only simulate flux during analysis period
    filedata = filedata.replace('mask_path = data/mask/<YYYY><MM>_maskfile.nc', 'mask_path ='+storm_dir_relative+ '/Resampled_for_fluxengine_storm_timings_with_masks.nc')
    filedata = filedata.replace('mask_prod = mask_variable','mask_prod = analysis_mask')
    filedata = filedata.replace('mask_temporalChunking = chunk_size','mask_temporalChunking ='+ time_chunk_val)

    #note that mask_timeDimensionName already set in config file template
    #mask_timeDimensionName = time

    # Write the file out again
    with open(configfilenew, 'w') as file:
      file.write(filedata)
    return configfilenew;

# if __name__ == "__main__":
def MAXSS_flux_run(MAXSS_working_directory="E:/MAXSS_working_directory",configfiletemplate="E:/MAXSS_Wilson/MAXSS_configuration_file_template.conf",verbose = True, specified_storms = [], test_run = False):
    # verbose=True
    #### Get the path of the root directory where the data are stored.
    #This will be user specific and can be changed depening on where data is stored
    # MAXSS_working_directory = "E:/MAXSS_working_directory";


    #note to use the same file structure used by the project r.g. #maxss/storm-atlas/ibtracts/region/year/storm
    os.chdir(MAXSS_working_directory);
    print("Working directory is now:", os.getcwd());

    #list of all the basins in MAXSS
    #MAXSS_regions=["east-pacific","north-atlantic","north-indian","south-atlantic","south-indian","south-pacific","west-pacific"]
    MAXSS_regions=["north-atlantic"]

#### Loop through the regions in MAXSS storm dataset
    for region in MAXSS_regions:

        #define the directory for the region
        region_directory = path.join(MAXSS_working_directory, "maxss", "storm-atlas", "tropical", "ibtracs", region);
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

            
            MAXSS_storms=storm_list
                
            storm_counter=0
            #### Loop through the storms for each year in the MAXSS storm dataset
            for storm in MAXSS_storms:

                # 1. If specified_storms has entries, check if ANY of your fragments match the storm name.
                if len(specified_storms) > 0:
                    # This checks if NONE of the fragments are inside the storm string
                    if not any(fragment in storm for fragment in specified_storms):
                        print(f"Skipping storm (no fragment match): {storm}")
                        storm_counter += 1  # Crucial: keep your tracking counter synchronized!
                        continue
                    
                print(f"--> Processing storm: {storm}")

                #directory for storm being processes
                storm_dir=storm_directory_list[storm_counter]
                #directory for storm being processes relative to current working directory
                storm_dir_relative = path.join("maxss", "storm-atlas", "tropical", "ibtracs", region, year, storm);

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
                with nc.Dataset(path.join("maxss", "storm-atlas", "tropical", "ibtracs", region, year, storm, "Resampled_for_fluxengine_MAXSS_L4_windspeed.nc")) as winds_nc:
                    wind_northward = winds_nc.variables['windspeed'][:]

                    wind_time_dimension=len(wind_northward)
                    timestepsinfile=wind_time_dimension

                    #### get start and end times - to be added to run call
                    wind_time = winds_nc.variables['time'][:]
                    wind_dates = num2date(wind_time, winds_nc.variables['time'].units)

                    time_chunk_val = str(len(wind_time))

                # Extract the land fraction
                with nc.Dataset(path.join("maxss", "storm-atlas", "tropical", "ibtracs", region, year, storm, "Resampled_for_fluxengine_MAXSS_land_fraction.nc")) as land_fraction_nc:
                    storm_land_fraction = land_fraction_nc.variables['land_proportion'][0]

                # Set model start and end times
                run_startime=wind_dates[0].strftime("%Y-%m-%d %H:%M")#
                
                # Check to see if this is a test run or full run
                if test_run:
                    # 5 days = 120 hourly steps. Cap at max available steps if storm is shorter than 5 days.
                    test_index = min(120, len(wind_dates) - 1)
                    run_endtime = wind_dates[test_index].strftime("%Y-%m-%d %H:%M")
                    print(f"[TEST RUN ACTIVE]: Restricting model timeframe to 5 days ({run_startime} to {run_endtime})")
                else:
                    run_endtime = wind_dates[-1].strftime("%Y-%m-%d %H:%M")

                #### Run flux engine for 'MAXSS run'
                run_name="MAXSS_RUN"
                # create custom config file for this storm
                # call custom function which copies file template and makes edits
                configFilePath_MAXSS_RUN=make_configuration_file(storm_dir_relative,timestepsinfile,region,year,storm,run_name,time_chunk_val,configfiletemplate)
                print("Running FluxEngine for Region={0} year={1} Storm={2}".format(region,year,storm));
                runStatus, fe_MAXSS_RUN = run_fluxengine(configFilePath_MAXSS_RUN,run_startime,run_endtime,processLayersOff=True, verbose=False);
                
                # Code to exit if test model run
                if test_run:
                    print(f"[TEST RUN COMPLETE]: Successfully verified 'MAXSS_RUN' for {storm}. Exiting storm loop as requested.")
                    return  # <exits the whole function

                #### Run flux engine for "REF run"
                run_name="REF_RUN"
                # create custom config file for this storm
                # call custom function which copies file template and makes edits
                configFilePath_REF_RUN=make_configuration_file(storm_dir_relative,timestepsinfile,region,year,storm,run_name,time_chunk_val,configfiletemplate)
                print("Running FluxEngine for Region={0} year={1} Storm={2}".format(region,year,storm));
                runStatus, fe_REF_RUN = run_fluxengine(configFilePath_REF_RUN,run_startime,run_endtime,processLayersOff=True, verbose=False);
                
                #### Run flux engine for "WIND run"
                run_name="WIND_RUN"
                # create custom config file for this storm
                # call custom function which copies file template and makes edits
                configFilePath_WIND_RUN=make_configuration_file(storm_dir_relative,timestepsinfile,region,year,storm,run_name,time_chunk_val,configfiletemplate)
                print("Running FluxEngine for Region={0} year={1} Storm={2}".format(region,year,storm));
                runStatus, fe_WIND_RUN = run_fluxengine(configFilePath_WIND_RUN,run_startime,run_endtime,processLayersOff=True, verbose=False);
                
                #### Run flux engine for "SST_with_gradients run"
                run_name="SST_WITH_GRADIENTS_RUN"
                # create custom config file for this storm
                # call custom function which copies file template and makes edits
                configFilePath_SST_WITH_GRADIENTS_RUN=make_configuration_file(storm_dir_relative,timestepsinfile,region,year,storm,run_name,time_chunk_val,configfiletemplate)
                print("Running FluxEngine for Region={0} year={1} Storm={2}".format(region,year,storm));
                runStatus, fe_SST_WITH_GRADIENTS_RUN = run_fluxengine(configFilePath_SST_WITH_GRADIENTS_RUN,run_startime,run_endtime,processLayersOff=True, verbose=False);
                
                #### Run flux engine for "SST_NO_gradients run"
                run_name="SST_NO_GRADIENTS_RUN"
                # create custom config file for this storm
                # call custom function which copies file template and makes edits
                configFilePath_SST_NO_GRADIENTS_RUN=make_configuration_file(storm_dir_relative,timestepsinfile,region,year,storm,run_name,time_chunk_val,configfiletemplate)
                print("Running FluxEngine for Region={0} year={1} Storm={2}".format(region,year,storm));
                runStatus, fe_SST_NO_GRADIENTS_RUN = run_fluxengine(configFilePath_SST_NO_GRADIENTS_RUN,run_startime,run_endtime,processLayersOff=True, verbose=False);
                
                #### Run flux engine for "SST_NO_gradients run"
                run_name="V_GAS_RUN"
                # create custom config file for this storm
                # call custom function which copies file template and makes edits
                configFilePath_V_GAS_RUN=make_configuration_file(storm_dir_relative,timestepsinfile,region,year,storm,run_name,time_chunk_val,configfiletemplate)
                print("Running FluxEngine for Region={0} year={1} Storm={2}".format(region,year,storm));
                runStatus, fe_V_GAS_RUN = run_fluxengine(configFilePath_V_GAS_RUN,run_startime,run_endtime,processLayersOff=True, verbose=False);
                
                #### Run flux engine for "SSS run"
                run_name="SSS_RUN"
                # create custom config file for this storm
                # call custom function which copies file template and makes edits
                configFilePath_SSS_RUN=make_configuration_file(storm_dir_relative,timestepsinfile,region,year,storm,run_name,time_chunk_val,configfiletemplate)
                print("Running FluxEngine for Region={0} year={1} Storm={2}".format(region,year,storm));
                runStatus, fe_SSS_RUN = run_fluxengine(configFilePath_SSS_RUN,run_startime,run_endtime,processLayersOff=True, verbose=False);
                
                #### Run flux engine for "PRESSURE run"
                run_name="PRESSURE_RUN"
                # create custom config file for this storm
                # call custom function which copies file template and makes edits
                configFilePath_PRESSURE_RUN=make_configuration_file(storm_dir_relative,timestepsinfile,region,year,storm,run_name,time_chunk_val,configfiletemplate)
                print("Running FluxEngine for Region={0} year={1} Storm={2}".format(region,year,storm));
                runStatus, fe_PRESSURE_RUN = run_fluxengine(configFilePath_PRESSURE_RUN,run_startime,run_endtime,processLayersOff=True, verbose=False);
                
                # #### Run flux engine for "PRECIPITATION run"
                # run_name="PRECIPITATION_RUN"
                # # create custom config file for this storm
                # # call custom function which copies file template and makes edits
                # configFilePath_PRECIPITATION_RUN=make_configuration_file(storm_dir_relative,timestepsinfile,region,year,storm,run_name)
                # print("Running FluxEngine for Region={0} year={1} Storm={2}".format(region,year,storm));
                # runStatus, fe_PRECIPITATION_RUN = run_fluxengine(configFilePath_PRECIPITATION_RUN,run_startime,run_endtime,processLayersOff=True, verbose=False);
                
                # Add to storm counter when everything is done.
                storm_counter=storm_counter+1

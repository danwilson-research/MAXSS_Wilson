#This script does all the plotting

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

runs=["MAXSS_RUN","REF_RUN","WIND_RUN","SST_RUN","SSS_RUN","PRESSURE_RUN"]
#runs=["MAXSS_RUN","REF_RUN","WIND_RUN","SST_RUN","PRESSURE_RUN"]

# 1. Define your color map at the top of your script or before the plot
run_colors = {
    "MAXSS_RUN": "darkgrey",
    "REF_RUN": "moccasin",
    "WIND_RUN": "k",
    "SST_RUN": "r",
    "SSS_RUN": "darkorchid",
    "PRESSURE_RUN": "orange"
}

def get_datetime(secondsSince1970):
    baseDate = datetime(1970,1,1);
    return baseDate+timedelta(seconds=secondsSince1970);

if __name__ == "__main__":
    verbose=True
#### Get the path of the root directory where the data are stored. 
    #This will be user specific and can be changed depening on where data is stored
    MAXSS_working_directory = "E:/MAXSS_working_directory"; 
    
    #Ensure output folders exist ---
    output_dir = os.path.join("output", "plots")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")
    
    #note to use the same file structure used by the project r.g. #maxss/storm-atlas/ibtracts/region/year/storm
    os.chdir(MAXSS_working_directory);
    print("Working directory is now:", os.getcwd());
    
    #list of all the basins in MAXSS 
    #MAXSS_regions=["east-pacific","north-atlantic","north-indian","south-atlantic","south-indian","south-pacific","west-pacific"]
    MAXSS_regions=["north-atlantic"]

#### Loop through the regions in MAXSS storm dataset
    for region in MAXSS_regions:
        
        #define the directory for the region
        region_directory = path.join(MAXSS_working_directory+r"\\maxss\\storm-atlas\\ibtracs\\{0}".format(region));
        
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
            
            #create empty matrix for bar data
            a=len(MAXSS_storms)
            b=len(runs)
            bar_plot_data = np.empty((a,b), dtype=float);
            
            
            storm_counter=0
        #### Loop through the storms for each year in the MAXSS storm dataset 
            for storm in MAXSS_storms:

                #directory for storm being processes
                storm_dir=storm_directory_list[storm_counter]
                #directory for storm being processes relative to current working directory
                storm_dir_relative = path.join(r"maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}".format(region,year,storm));
        
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
                
                Timeseries_plot_data= [];
                runcounter=0
                for run_type in runs:
                    name_of_run=runs[runcounter]
                    # path to integrated fluxes 
                    int_flux_paths_nc=nc.Dataset(path.join(r"output\Spatially_integrated_fluxes\maxss\\storm-atlas\\ibtracs\\{0}\\{1}\\{2}_{3}.nc".format(region,year,storm,run_type)));
    
                    Hourlyflux_for_this_run= int_flux_paths_nc.variables['Hourly_flux'][:]
                    Hourlyflux_for_this_run[Hourlyflux_for_this_run==0] = np.nan
                    Hourlyfluxdate_for_this_run= int_flux_paths_nc.variables['time'][:]
                    Total_flux_for_this_run= int_flux_paths_nc.variables['Total_flux'][:]
                    


                    dates = [get_datetime(int(timeSecs)) for timeSecs in int_flux_paths_nc.variables["time"][:]];


                    #save data to 2d array to plot as bars
                    bar_plot_data[storm_counter,runcounter]=Total_flux_for_this_run
                    Timeseries_plot_data.append(Hourlyflux_for_this_run)

                
                
                    #### make timeseries figure of the flux for this storm
                    fig1 = plt.figure(figsize=(24,15))
                    gs = fig1.add_gridspec(1, 1)
                    f1_ax1 = fig1.add_subplot(gs[0, :])
                    plt.plot(dates,Hourlyflux_for_this_run, linewidth=1.5)
                    plt.ylabel("Total flux in region (Tg C hr${^-1}$)",fontsize=36)
                    plt.xlabel("Time",fontsize=36)
                    f1_ax1.xaxis.set_major_formatter(Month_Fmt)
                    plt.xticks(fontsize=36)
                    plt.yticks(fontsize=36)
                    plt.legend([name_of_run],fontsize=36)
                    f1_ax1.tick_params(labelcolor='k', direction='in',length=10,width=4)
                    plt.show()

                    runcounter=runcounter+1


                # all of the plots where we want to compare runs
                fig = plt.figure(figsize=(24,15))
                gs = fig.add_gridspec(1, 1)
                f_ax1 = fig.add_subplot(gs[0, :])
                barlist=f_ax1.bar(runs, bar_plot_data[0], color = 'b', width = 0.95)
                
                
                # 2. Create the bar plot
                barlist = f_ax1.bar(runs, bar_plot_data[storm_counter], color='b', width=0.95)
                
                # 3. Use a loop to color only the bars that exist
                for i, run_name in enumerate(runs):
                    # This checks if the run name has a specific color assigned
                    if run_name in run_colors:
                        barlist[i].set_color(run_colors[run_name])
                #plt.title("Total flux over storm domain for storm period",fontsize=36)
                plt.ylabel('Total flux (Tg C)',fontsize=36)
                plt.xlabel('Runs',fontsize=36)
                f_ax1.set_xticklabels(runs,rotation=90,fontsize=36)
                plt.yticks(fontsize=36)
                plt.show()
                fig.savefig(r"output\plots\\\Integrated_fluxes_for_storm_{0}.png".format(storm));
    
    
                    # all of the plots where we want to compare runs
                fig = plt.figure(figsize=(24,15))
                gs = fig.add_gridspec(1, 1)
                f_ax1 = fig.add_subplot(gs[0, :])
                barlist=f_ax1.bar(["MAXSS_RUN","REF_RUN"], bar_plot_data[0][0:2], color = 'b', width = 0.95)
                barlist[0].set_color('darkgrey')#main
                barlist[1].set_color('moccasin')#ref
                # barlist[2].set_color('k')#wind
                # barlist[3].set_color('r')#sst
                # barlist[4].set_color('darkorchid')#sss
                # barlist[5].set_color('orange')#pressure
                #plt.title("Total flux over storm domain for storm period",fontsize=36)
                plt.ylabel('Total flux (Tg C)',fontsize=36)
                plt.xlabel('Runs',fontsize=36)
                f_ax1.set_xticklabels(runs,rotation=90,fontsize=36)
                plt.yticks(fontsize=36)
                plt.show()
                fig.savefig(r"output\plots\\\Integrated_fluxes_for_storm_{0}_mainandref.png".format(storm));
    
    
                fig1 = plt.figure(figsize=(24,15))
                gs = fig1.add_gridspec(1, 1)
                f1_ax1 = fig1.add_subplot(gs[0, :])
                plt.plot(dates,Timeseries_plot_data[0],color="gray", linewidth=3)
                plt.ylabel("Total flux in region (Tg C hr${^-1}$)",fontsize=36)
                plt.xlabel("Time",fontsize=36)
                f1_ax1.xaxis.set_major_formatter(Month_Fmt)
                plt.xticks(fontsize=36)
                plt.yticks(fontsize=36)
                plt.legend(runs,fontsize=36)
                plt.plot(dates,np.zeros(len(dates)),"--",color="k", linewidth=3)
                f1_ax1.tick_params(labelcolor='k', direction='in',length=10,width=4)
                plt.show()
                fig1.savefig(r"output\plots\\\Flux_timeseries_for_storm_{0})main_run.png".format(storm));

                
                fig1 = plt.figure(figsize=(24,15))
                gs = fig1.add_gridspec(1, 1)
                f1_ax1 = fig1.add_subplot(gs[0, :])
                plt.plot(dates,Timeseries_plot_data[0],color="gray", linewidth=3)
                plt.plot(dates,Timeseries_plot_data[1],color="brown", linewidth=3)
                plt.legend(runs,fontsize=36)
                plt.plot(dates,np.zeros(len(dates)),"--",color="k", linewidth=3)
                plt.ylabel("Total flux in region (Tg C hr${^-1}$)",fontsize=36)
                plt.xlabel("Time",fontsize=36)
                f1_ax1.xaxis.set_major_formatter(Month_Fmt)
                plt.xticks(fontsize=36)
                plt.yticks(fontsize=36)
                f1_ax1.tick_params(labelcolor='k', direction='in',length=10,width=4)
                plt.show()
                fig1.savefig(r"output\plots\\\Flux_timeseries_for_storm_{0}main_run_and_ref.png".format(storm));

                
    
                fig1 = plt.figure(figsize=(24,15))
                gs = fig1.add_gridspec(1, 1)
                f1_ax1 = fig1.add_subplot(gs[0, :])
                plt.plot(dates,Timeseries_plot_data[0],color="gray", linewidth=3)
                plt.plot(dates,Timeseries_plot_data[1],color="brown", linewidth=3)
                plt.plot(dates,Timeseries_plot_data[2],color="k", linewidth=3)
                plt.plot(dates,Timeseries_plot_data[3],color="r", linewidth=3)
                plt.plot(dates,Timeseries_plot_data[4],color="darkorchid", linewidth=3)
                plt.plot(dates,Timeseries_plot_data[5],color="orange", linewidth=3)
                plt.plot(dates,np.zeros(len(dates)),"--",color="k", linewidth=3)
                plt.ylabel("Total flux in region (Tg C hr${^-1}$)",fontsize=36)
                plt.xlabel("Time",fontsize=36)
                f1_ax1.xaxis.set_major_formatter(Month_Fmt)
                plt.xticks(fontsize=36)
                plt.yticks(fontsize=36)
                plt.legend(runs,fontsize=36)
                f1_ax1.tick_params(labelcolor='k', direction='in',length=10,width=4)
                plt.show()
                fig1.savefig(r"output\plots\\\Flux_timeseries_for_storm_{0}.png".format(storm));

                
                
                # Add to storm counter when everything is done.
                storm_counter=storm_counter+1
                
            #now do any plots of all the storms in a region
            # Bar graph of the integrated flux by storm 
            
            
            
            # fig = plt.figure()
            # ax = fig.add_axes([0,0,1,1])
            # barlist=ax.bar(runs, bar_plot_data[0], color = 'b', width = 0.95)
            # barlist[0].set_color('darkgrey')#main
            # barlist[1].set_color('moccasin')#ref
            # barlist[2].set_color('k')#wind
            # barlist[3].set_color('r')#sst
            # barlist[4].set_color('darkorchid')#sss
            # barlist[5].set_color('orange')#pressure
            # plt.title("Total flux over storm domain for storm period")
            # plt.ylabel('Total flux (Tg C)')
            # plt.xlabel('Runs')
            # plt.show()








                

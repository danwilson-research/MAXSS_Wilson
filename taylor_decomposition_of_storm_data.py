# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 11:25:01 2026

@author: dw557
"""

#Script to run a taylor decomposition for each storm.

#Import required packages
import os
from os import path
from glob import glob
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

#List out model components
total_flux_run = "MAXSS_RUN"
reference_run = "REF_RUN"
wind_run = "WIND_RUN"
sst_run = "SST_RUN"
sss_run ="SSS_RUN"
pressure_run = "PRESSURE_RUN"
precipitation_run ="PRECIPITATION_RUN"

MAXSS_working_directory = "E:/MAXSS_working_directory/"

plot_save_location = "E:/MAXSS_working_directory/taylor_decomposition/"

#CHANGE WHEN I MOVE TO INLCUDING MORE REGIONS
MAXSS_regions = ["north-atlantic"]

# --- MAIN LOOP ---
for region in MAXSS_regions: 
    
    # Define the directory for the region
    region_directory = path.join(MAXSS_working_directory, "output/spatially_integrated_fluxes/maxss/storm-atlas/ibtracs", region)
    
    # Get all year directories
    year_directory_list = glob(path.join(region_directory, "*/"))
    
    for year_dir in year_directory_list:
        # Extract year name from path for plot labelling
        year_name = os.path.basename(os.path.normpath(year_dir))
        
        # Get a list of all 'Total' runs to identify unique storms
        total_run_files = glob(os.path.join(year_dir, f"*{total_flux_run}.nc"))
        
        for total_file in total_run_files:
            #Extract the storm prefix (e.g., 010215N13319_AL042010_COLIN)
            file_base = os.path.basename(total_file)
            storm_id = file_base.replace(f"_{total_flux_run}.nc", "")
       
            # Load the .nc files using xarray
            # Use the unique storm_id to target the matching drivers
            try:
                ds_total = xr.open_dataset(total_file)
                ds_ref  = xr.open_dataset(os.path.join(year_dir, f"{storm_id}_{reference_run}.nc"))
                ds_wind = xr.open_dataset(os.path.join(year_dir, f"{storm_id}_{wind_run}.nc"))
                ds_sst  = xr.open_dataset(os.path.join(year_dir, f"{storm_id}_{sst_run}.nc"))
                ds_sss  = xr.open_dataset(os.path.join(year_dir, f"{storm_id}_{sss_run}.nc"))
                ds_pressure  = xr.open_dataset(os.path.join(year_dir, f"{storm_id}_{pressure_run}.nc"))
                ds_precipitation  = xr.open_dataset(os.path.join(year_dir, f"{storm_id}_{precipitation_run}.nc"))
            
                # Calculate Individual Contributions (Taylor Terms)
                # F_driver - F_ref
                wind_contrib   = ds_wind['Hourly_flux'] - ds_ref['Hourly_flux']
                sst_contrib    = ds_sst['Hourly_flux'] - ds_ref['Hourly_flux']
                sss_contrib    = ds_sss['Hourly_flux'] - ds_ref['Hourly_flux']
                pres_contrib   = ds_pressure['Hourly_flux'] - ds_ref['Hourly_flux']
                precip_contrib = ds_precipitation['Hourly_flux'] - ds_ref['Hourly_flux']
                
                # Sum the Taylor components 
                taylor_sum = wind_contrib + sst_contrib + sss_contrib + pres_contrib + precip_contrib
                
                # Calculate the actual retrieved anomaly (The "Truth")
                total_anomaly = ds_total['Hourly_flux'] - ds_ref['Hourly_flux']

                #calculate residual between total anomaly and Taylor sum
                residual = total_anomaly - taylor_sum
                
                ## Plot Taylor decomposition breakdown ##
                
                # Convert xarray DataArrays to pandas for easier plotting
                time_axis = ds_total.time.values
                
                plt.figure(figsize=(14, 8))
                
                # Plot Individual Taylor Components (The Drivers)
                plt.plot(time_axis, wind_contrib, label='Wind Speed Contribution', color='black', linewidth=1.5)
                plt.plot(time_axis, sst_contrib, label='SST Contribution', color='red', linewidth=1.5)
                plt.plot(time_axis, sss_contrib, label='SSS Contribution', color='purple', alpha=0.6)
                plt.plot(time_axis, precip_contrib, label='Precipitation Contribution', color='blue', alpha=0.6)
                plt.plot(time_axis, pres_contrib, label='Pressure Contribution', color='orange', alpha=0.6)
                
                # Plot the Taylor Sum vs. The Actual Total Anomaly
                plt.plot(time_axis, taylor_sum, label='Taylor Sum ', 
                         color='cyan', linestyle='--', linewidth=2, alpha=0.8)
                plt.plot(time_axis, total_anomaly, label='Actual Total Anomaly (MAXSS_RUN)', 
                         color='grey', linestyle=':', linewidth=2)
                
                #  Formatting the Plot
                plt.axhline(0, color='black', lw=1, ls='-') # Zero baseline
                plt.title(f"Taylor Decomposition of Sea-Air CO2 Flux Anomalies\nStorm: {storm_id} ({year_name})")
                plt.xlabel("Time")
                plt.ylabel("Flux Anomaly (Tg C hr⁻¹)")
                plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
                plt.grid(True, which='both', linestyle='--', alpha=0.5)
                plt.tight_layout()
                
                # Save or Show
                plt.savefig(path.join(plot_save_location, f"hrly_taylor_decomp_{storm_id}.png"), dpi=300)
                plt.show()
                
                ## Plot Taylor decomposition breakdown contribution to sea-air flux over whoel study period ##
                
                # Sum the hourly rates to get total Tg C for the 60 days
                # data is hourly, so the sum of Tg C/hr is the total Tg C
                driver_totals = {
                    'Wind': wind_contrib.sum().values,
                    'SST': sst_contrib.sum().values,
                    'SSS': sss_contrib.sum().values,
                    'Precip': precip_contrib.sum().values,
                    'Pressure': pres_contrib.sum().values,
                    'Taylor Sum': taylor_sum.sum().values,
                    'Actual Total': total_anomaly.sum().values }
                                
                # Create the Bar Plot
                labels = list(driver_totals.keys())
                values = list(driver_totals.values())
                colors = ['black', 'firebrick', 'rebeccapurple', 'royalblue', 'goldenrod', 'cyan', 'grey']
                
                plt.figure(figsize=(10, 6))
                bars = plt.bar(labels, values, color=colors, alpha=0.8)
                
                # Add labels and styling
                plt.axhline(0, color='black', lw=1)
                plt.ylabel('Total Carbon Flux Anomaly (Tg C)')
                plt.title(f'Cumulative Driver Contributions \nStorm: {storm_id}')
                plt.xticks(rotation=45)
                
                # Add value labels on top of bars
                for bar in bars:
                    yval = bar.get_height()
                    plt.text(bar.get_x() + bar.get_width()/2, yval, round(float(yval), 4), 
                             va='bottom' if yval > 0 else 'top', ha='center', fontsize=9)
                
                plt.tight_layout()
                plt.savefig(path.join(plot_save_location, f"cumulative_impact_{storm_id}.png"), dpi=300)
                plt.show()
                
            except FileNotFoundError:
                print(f"Missing a component for {storm_id}, skipping...")
                continue
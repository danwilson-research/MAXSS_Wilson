# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 11:25:01 2026

@author: dw557
"""

#Script to run a taylor decomposition for each storm.

#Import required packages
# Script to run a taylor decomposition for each storm.

import os
from os import path
from glob import glob
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

# Model Run Definitions
total_flux_run = "MAXSS_RUN"
reference_run = "REF_RUN"
wind_run = "WIND_RUN"
sst_no_grad_run = "SST_NO_GRADIENTS_RUN"
sst_with_grad_run = "SST_WITH_GRADIENTS_RUN"
sss_run = "SSS_RUN"
v_gas_run = "V_GAS_RUN"
pressure_run = "PRESSURE_RUN"

MAXSS_working_directory = "E:/MAXSS_working_directory/"
plot_save_location = "E:/MAXSS_working_directory/output/plots/taylor_decomposition/"
os.makedirs(plot_save_location, exist_ok=True)

MAXSS_regions = ["north-atlantic"]
all_storm_data = []

# Main code loop
for region in MAXSS_regions: 
    region_directory = path.join(MAXSS_working_directory, "output/spatially_integrated_fluxes/maxss/storm-atlas/ibtracs", region)
    year_directory_list = glob(path.join(region_directory, "*/"))
    
    for year_dir in year_directory_list:
        year_name = os.path.basename(os.path.normpath(year_dir))
        total_run_files = glob(os.path.join(year_dir, f"*{total_flux_run}.nc"))
        
        for total_file in total_run_files:
            file_base = os.path.basename(total_file)
            storm_id = file_base.replace(f"_{total_flux_run}.nc", "")
       
            try:
                # Load Datasets
                ds_total = xr.open_dataset(total_file)
                ds_ref   = xr.open_dataset(os.path.join(year_dir, f"{storm_id}_{reference_run}.nc"))
                ds_wind  = xr.open_dataset(os.path.join(year_dir, f"{storm_id}_{wind_run}.nc"))
                ds_sst_no_grad   = xr.open_dataset(os.path.join(year_dir, f"{storm_id}_{sst_no_grad_run}.nc"))
                ds_sst_with_grad = xr.open_dataset(os.path.join(year_dir, f"{storm_id}_{sst_with_grad_run}.nc"))
                ds_sss   = xr.open_dataset(os.path.join(year_dir, f"{storm_id}_{sss_run}.nc"))
                ds_v_gas = xr.open_dataset(os.path.join(year_dir, f"{storm_id}_{v_gas_run}.nc"))
                ds_pressure = xr.open_dataset(os.path.join(year_dir, f"{storm_id}_{pressure_run}.nc"))
            
                # Calculate Base Contributions
                wind_c  = ds_wind['Hourly_flux'] - ds_ref['Hourly_flux']
                sss_c   = ds_sss['Hourly_flux'] - ds_ref['Hourly_flux']
                vgas_c  = ds_v_gas['Hourly_flux'] - ds_ref['Hourly_flux']
                pres_c  = ds_pressure['Hourly_flux'] - ds_ref['Hourly_flux']
                
                # Methodology A: No Gradients
                sst_no_g_c = ds_sst_no_grad['Hourly_flux'] - ds_ref['Hourly_flux']
                taylor_sum_no_g = wind_c + sst_no_g_c + sss_c + vgas_c + pres_c
                
                # Methodology B: With Gradients
                sst_with_g_c = ds_sst_with_grad['Hourly_flux'] - ds_ref['Hourly_flux']
                taylor_sum_with_g = wind_c + sst_with_g_c + sss_c + vgas_c + pres_c
                
                total_anomaly = ds_total['Hourly_flux'] - ds_ref['Hourly_flux']
                
                # Masking and Trimming
                valid_mask = (total_anomaly != 0)
                last_idx = np.where(valid_mask.values)[0][-1] if any(valid_mask) else None
                
                if last_idx is None:
                    continue

                # Prepare Plotting Slices
                t_plot = ds_total.time.values[:last_idx]
                
                # Define Scenarios for looping
                scenarios = [
                    {'name': 'SST No Gradients', 'sst_c': sst_no_g_c, 'sum': taylor_sum_no_g, 'color': 'blue'},
                    {'name': 'SST With Gradients', 'sst_c': sst_with_g_c, 'sum': taylor_sum_with_g, 'color': 'cyan'}
                ]

                for sc in scenarios:
                    # 1. Hourly Timeseries Plot
                    plt.figure(figsize=(12, 6))
                    plt.plot(t_plot, wind_c[:last_idx], label='Wind', color='grey', alpha=0.7)
                    plt.plot(t_plot, sc['sst_c'][:last_idx], label=f'SST ({sc["name"]})', color=sc['color'], linewidth=2)
                    plt.plot(t_plot, sss_c[:last_idx], label='SSS', color='purple', alpha=0.5)
                    plt.plot(t_plot, vgas_c[:last_idx], label='V Gas', color='green', alpha=0.5)
                    plt.plot(t_plot, pres_c[:last_idx], label='Pressure', color='saddlebrown', alpha=0.5)
                    
                    plt.plot(t_plot, sc['sum'][:last_idx], label='Taylor Sum', color='magenta', linestyle='--')
                    plt.plot(t_plot, total_anomaly[:last_idx], label='Actual Anomaly', color='black', linestyle=':', linewidth=2)
                    
                    plt.title(f"Taylor Decomposition ({sc['name']})\nStorm: {storm_id}")
                    plt.ylabel("Flux Anomaly (Tg C hr⁻¹)")
                    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                    plt.grid(True, alpha=0.3)
                    plt.tight_layout()
                    plt.savefig(path.join(plot_save_location, f"timeseries_{sc['name']}_{storm_id}.png"))
                    plt.close()

                    # 2. Cumulative Bar Plot
                    drv_labels = ['Wind', 'SST', 'SSS', 'V Gas', 'Pres', 'Taylor Sum', 'Actual Anomaly']
                    drv_values = [
                        wind_c[:last_idx].sum().item(),
                        sc['sst_c'][:last_idx].sum().item(),
                        sss_c[:last_idx].sum().item(),
                        vgas_c[:last_idx].sum().item(),
                        pres_c[:last_idx].sum().item(),
                        sc['sum'][:last_idx].sum().item(),
                        total_anomaly[:last_idx].sum().item()]
                    
                    fig, ax = plt.subplots(figsize=(10, 6))
                    colors = ['grey', sc['color'], 'purple', 'green', 'saddlebrown', 'magenta', 'black']
                    bars = ax.bar(drv_labels, drv_values, color=colors, alpha=0.8)
                    
                    # Add labels to the top of each bar ---
                    ax.bar_label(bars, fmt='%.3f', padding=3, fontsize=9)
                    
                    # Styling
                    ax.axhline(0, color='black', lw=0.8)
                    ax.set_title(f"Cumulative Impact ({sc['name']})\nStorm: {storm_id}")
                    ax.set_ylabel("Total Flux Anomaly (Tg C)")
                    
                    # Increase y-limit slightly so labels don't get cut off at the top
                    ax.set_ylim(ax.get_ylim()[0] * 1.1, ax.get_ylim()[1] * 1.1)
                    
                    plt.xticks(rotation=45, ha='right')
                    plt.tight_layout()
                    plt.savefig(path.join(plot_save_location, f"bar_{sc['name']}_{storm_id}.png"), dpi=300)
                    plt.close()

                # Save Data to CSV list
                # We calculate the sums once here to keep the dictionary clean
                all_storm_data.append({
                    'Storm_ID': storm_id, 
                    'Year': year_name,
                    'Region': region,
                    # Individual Drivers (Common to both methods)
                    'Wind_Contribution_TgC': wind_c[:last_idx].sum().item(),
                    'SSS_Contribution_TgC': sss_c[:last_idx].sum().item(),
                    'VGas_Contribution_TgC': vgas_c[:last_idx].sum().item(),
                    'Pressure_Contribution_TgC': pres_c[:last_idx].sum().item(),
                    # SST Specific Drivers
                    'SST_NoGrad_Contribution_TgC': sst_no_g_c[:last_idx].sum().item(),
                    'SST_WithGrad_Contribution_TgC': sst_with_g_c[:last_idx].sum().item(),
                    # Totals and Verification
                    'Actual_Anomaly_TgC': total_anomaly[:last_idx].sum().item(),
                    'Taylor_Sum_NoGrad_TgC': taylor_sum_no_g[:last_idx].sum().item(),
                    'Taylor_Sum_WithGrad_TgC': taylor_sum_with_g[:last_idx].sum().item()
                })
                
            except FileNotFoundError:
                continue

# Save Summary CSV
if all_storm_data:
    # 1. Define the directory and filename
    # Removed the leading slash to ensure it stays inside MAXSS_working_directory
    csv_dir = path.join(MAXSS_working_directory, "output", "taylor_decomposition_summary")
    csv_file = path.join(csv_dir, "taylor_decomposition_summary.csv")
    
    # 2. Create the directory if it doesn't exist
    os.makedirs(csv_dir, exist_ok=True)
    
    # 3. Save the file
    df_summary = pd.DataFrame(all_storm_data)
    df_summary.to_csv(csv_file, index=False)
    
    print(f"CSV summary saved to: {csv_file}")
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 16:15:25 2026

@author: dw557
"""

# Script to compare air-sea co2 flux between Dan Wilson FluxENgine output and 
# 'UExP-FNN-U full surface ocean carbonate system' data
# Updated to account for updated Fluxengine output that tracks the storm 

# --- SCRIPT: AIR-SEA CO2 FLUX COMPARISON ---
# Purpose: Compare FluxEngine output (tracking storms) with UExP global datasets.
# Logic: Hourly regridding and common masking to ensure robust comparison.

import xarray as xr
import numpy as np
import os
import pandas as pd
from os import path
from glob import glob
import matplotlib.pyplot as plt

#NEED TO ACCOUNT FOR LAND PROPORTION IN FLUXENGINE AND DAN F DATA ##

# --- CONFIGURATION ---
MAXSS_working_directory = "E:/MAXSS_working_directory"
output_base = 'E:/MAXSS_working_directory/FluxEngine_UExP_comparison/'
individual_storm_dir = os.path.join(output_base, 'individual_storm_data')
UEXP_FNN_path = "E:/MAXSS_working_directory/Ford_et_al_GBC_fco2/Fordetal_UExP-FNN-U_surface-carbonate-system_v2025-1.nc"

# Ensure output directories exist
for d in [output_base, individual_storm_dir]:
    if not os.path.exists(d):
        os.makedirs(d)

# Parameters for regions and skip list
MAXSS_regions = ["north-atlantic"]
storms_to_skip = [] #"BONNIE", "COLIN", "MARIA", "RINA"

# Load Global UExP Data (Opened once to save memory)
UExP_FNN_raw = xr.open_dataset(UEXP_FNN_path)

# List to hold summary data (one row per storm)
summary_data = []

# --- MAIN LOOP ---
for region in MAXSS_regions:
    region_directory = path.join(MAXSS_working_directory, "output/MAXSS_RUN/maxss/storm-atlas/ibtracs", region)
    year_directory_list = glob(path.join(region_directory, "*/"))
    
    for year_dir in year_directory_list:
        year_name = os.path.basename(os.path.normpath(year_dir))
        storm_directory_list = glob(path.join(year_dir, "*/"))

        for storm_path in storm_directory_list:
            storm_name = os.path.basename(os.path.normpath(storm_path))
            
            if any(name in storm_name for name in storms_to_skip):
                print(f"Skipping storm: {storm_name}")
                continue
            
            print(f'Processing: {storm_name} ({year_name})')
            
            fluxengine_dir = storm_path
            fluxengine_files = [f for f in os.listdir(fluxengine_dir) if f.endswith('.nc')]
            
            if not fluxengine_files:
                print(f"Skipping {storm_name}: No NetCDF files found.")
                continue

            # List to hold hourly data for THIS SPECIFIC STORM
            current_storm_hourly_data = []

            for file in fluxengine_files:
                file_path = os.path.join(fluxengine_dir, file)
                
                with xr.open_dataset(file_path) as fluxengine_hourly_ds:
                    for t in range(len(fluxengine_hourly_ds.time)):
                        # 1. Slice and Time extraction
                        fe_flux_slice = fluxengine_hourly_ds.OF.isel(time=t)
                        current_time = fe_flux_slice.time.values
                    
                        # 2. Get matching UExP slice
                        uexp_slice = UExP_FNN_raw.flux.sel(time=current_time, method="nearest").transpose("latitude", "longitude")
            
                        # 3. Regrid UExP to FluxEngine grid (0.25 deg)
                        uexp_regridded = uexp_slice.interp_like(fe_flux_slice, method="nearest")
            
                        # 4. Create common mask (intersection of both datasets)
                        common_mask = fe_flux_slice.notnull() & uexp_regridded.notnull()
            
                        # 5. Apply mask
                        fe_final = fe_flux_slice.where(common_mask)
                        uexp_final = uexp_regridded.where(common_mask)
                        
                        # 6. Integrated Hourly Flux Calculation (Tg C hr-1)
                        R = 6371000  # Earth's radius
                        dlat_deg = 0.25
                        dlon_rad = np.radians(0.25)
                        
                        lat_top_rad = np.radians(fe_final.latitude + (dlat_deg / 2))
                        lat_bot_rad = np.radians(fe_final.latitude - (dlat_deg / 2))
                        
                        pixel_areas = (R**2) * dlon_rad * (np.sin(lat_top_rad) - np.sin(lat_bot_rad))

                        # Summing Mass = (Rate * Area) / 24 hrs
                        fe_sum_g_day = (fe_final * pixel_areas).sum(dim=['latitude', 'longitude'])
                        uexp_sum_g_day = (uexp_final * pixel_areas).sum(dim=['latitude', 'longitude'])
                        
                        fe_total_hr_tg = (fe_sum_g_day * 1e-12) / 24.0
                        uexp_total_hr_tg = (uexp_sum_g_day * 1e-12) / 24.0

                        # 7. Area Stats for robustness tracking
                        comparison_area_km2 = (common_mask * pixel_areas).sum() / 1e6
                        area_lost_mask = fe_flux_slice.notnull() & ~common_mask
                        area_lost_km2 = (area_lost_mask * pixel_areas).sum() / 1e6

                        # Store hourly slice result
                        current_storm_hourly_data.append({
                            'time': str(current_time),
                            'fe_tg_hr': float(fe_total_hr_tg),
                            'uexp_tg_hr': float(uexp_total_hr_tg),
                            'comparison_area_km2': float(comparison_area_km2),
                            'area_lost_km2': float(area_lost_km2)})
                        
                        # Quick check print using scientific notation (Tg can be very small hourly)
                        #print(f"Time: {current_time} | FE: {fe_total_hr_tg.values:.4f} Tg/hr | UExP: {uexp_total_hr_tg.values:.4f} Tg/hr")
                        
                        
            # --- END OF STORM: SAVE HOURLY DATA ---
            if current_storm_hourly_data:
                df_hourly = pd.DataFrame(current_storm_hourly_data)
                hourly_filename = f"{storm_name}_{year_name}_hourly_flux.csv"
                df_hourly.to_csv(os.path.join(individual_storm_dir, hourly_filename), index=False)

                # Append summary row for the global summary file
                summary_data.append({
                    'region': region,
                    'storm_name': storm_name,
                    'year': year_name,
                    'total_fe_tg': df_hourly['fe_tg_hr'].sum(),
                    'total_uexp_tg': df_hourly['uexp_tg_hr'].sum(),
                    'mean_comparison_area_km2': df_hourly['comparison_area_km2'].mean(),
                    'total_area_lost_km2': df_hourly['area_lost_km2'].sum(),
                    'hours_processed': len(df_hourly)
                })

# --- END OF ALL LOOPS: SAVE GLOBAL SUMMARY ---
if summary_data:
    df_summary = pd.DataFrame(summary_data)
    df_summary.to_csv(os.path.join(output_base, 'all_storms_summary.csv'), index=False)
    print(f"Workflow Complete. Summary saved to: {output_base}")
else:
    print("No data was processed.")
                        
             
#Plot comparison bar chart per storm between UEXP and FluxEngine               
             
# 1. Extract just the storm name (the part after the last underscore)
# This creates a list of strings like ['ALEX', 'BONNIE', 'COLIN', ...]
clean_names = [name.split('_')[-1] for name in df['storm_name']]

# 2. Setup the data positions
x = np.arange(len(clean_names)) 
width = 0.35 

fig, ax = plt.subplots(figsize=(10, 6))

# 3. Create the grouped bars
rects1 = ax.bar(x - width/2, df['total_fe_tg'], width, label='FluxEngine', color='#3498db')
rects2 = ax.bar(x + width/2, df['total_uexp_tg'], width, label='UEXP-FNN', color='#e74c3c')

# 4. Customizing the labels
ax.set_ylabel('Total Gas (tg)')
ax.set_title('Comparison of Total FluxEngine vs UEXP ')

# Use the cleaned names for the x-axis ticks
ax.set_xticks(x)
ax.set_xticklabels(clean_names) 

ax.legend()
ax.grid(axis='y', linestyle='--', alpha=0.6)

fig.tight_layout()

plt.show()          
                
             
                
             
                


ALSO DOUBLE CHECK I KNOW WHAT MISSING AREA IS SHOWING AND COMPARE TO PREVIOUS PLOTS


#Following block of code used to sanity check 

test = uexp_final.data
print('uexp' + str(np.nansum(test)))

test1 = fe_final.data
print('FE' + str(np.nansum(test1)))
            
                
                
            
            
            
            
            
            
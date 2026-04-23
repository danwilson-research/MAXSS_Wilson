# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 17:30:54 2026

@author: dw557
"""

# Script to calculate the hourly fluxes throughout a storm

#Import required packages
import xarray as xr
import numpy as np
import os
import pandas as pd
from os import path
from glob import glob
from pyproj import Geod
from netCDF4 import Dataset
from tqdm import tqdm

# --- CONFIGURATION ---
MAXSS_working_directory = "E:/MAXSS_working_directory"
output_base = 'E:/MAXSS_working_directory/output/Spatially_integrated_fluxes'
netcdf_output_root = path.join(MAXSS_working_directory, "output/Spatially_integrated_fluxes/maxss/storm-atlas/ibtracs")


# The list of runs to process
runs = ["MAXSS_RUN", "REF_RUN", "WIND_RUN", "SST_NO_GRADIENTS_RUN", 
        "SST_WITH_GRADIENTS_RUN", "SSS_RUN", "V_GAS_RUN", "PRESSURE_RUN"]

# Regions and storms to process
MAXSS_regions = ["north-atlantic"]
storms_to_skip = ["ALEX"]  #, "BONNIE", "COLIN", "MARIA", "RINA"

# Use a dictionary to store data before converting to a wide CSV
# Format: {(Storm, Year, Region): {Run1: Total, Run2: Total}}
master_summary = {}

# Initialize the Geod object with the WGS84 ellipsoid
geod = Geod(ellps="WGS84")

def calculate_ellipsoidal_areas(ds, res=0.25):
    storm_lat = ds.latitude.values
    storm_lon = ds.longitude.values
    
    # We only need to calculate the area for one longitude to get the row value
    lon_sample = storm_lon[0]
    row_areas = []

    for lat in storm_lat:
        # Define the corners just like your working script
        lats = [lat-res/2, lat-res/2, lat+res/2, lat+res/2, lat-res/2]
        lons = [lon_sample-res/2, lon_sample+res/2, lon_sample+res/2, lon_sample-res/2, lon_sample-res/2]
        
        # Use the specific method from your working script
        poly_area, _ = geod.polygon_area_perimeter(lons, lats)
        row_areas.append(abs(poly_area))
        
    # Return as xarray for automatic broadcasting
    return xr.DataArray(row_areas, coords={'latitude': storm_lat}, dims='latitude')

# --- MAIN LOOP ---

# Loop through regions
for region in MAXSS_regions:
    # We use the MAXSS_RUN directory to establish the list of available storms
    reference_directory = path.join(MAXSS_working_directory, "output/MAXSS_RUN/maxss/storm-atlas/ibtracs", region)
    year_directory_list = glob(path.join(reference_directory, "*/"))
    
    # Loop through years of data
    for year_dir in year_directory_list:
        year_name = os.path.basename(os.path.normpath(year_dir))
        storm_directory_list = glob(path.join(year_dir, "*/"))

         
        # Loop through storms
        for storm_path in storm_directory_list:
            storm_name = os.path.basename(os.path.normpath(storm_path))
            
            # skip the following storms (helpful for debugging)
            if any(name in storm_name for name in storms_to_skip):
                continue
            
            land_fraction_path = path.join(MAXSS_working_directory, "maxss/storm-atlas/ibtracs/", region, year_name,storm_name,'Resampled_for_fluxengine_MAXSS_land_fraction.nc')
            
            #Load the land proportion mask
            with xr.open_dataset(land_fraction_path) as lf_ds:
                # Use .load() if the file is small to speed up indexing later
                land_fraction = lf_ds.land_proportion.load() 

            # Handle the FillValue (9.969e+36)
            # Any value larger than 1.0 is clearly a mask/fill value. 
            land_fraction = land_fraction.where(land_fraction <= 1.0, 1.0)

            #Extract the water proportion
            water_proportion = 1.0 - land_fraction
            
            # Hold areas for the storm, calculated on the first componnent run found
            grid_areas = None
            
            run_pbar = tqdm(runs, desc="Component Runs", leave=False)
            
            #Loop through component runs
            for run_name in run_pbar:
                run_pbar.set_description(f"{storm_name}: {run_name}")
            
                # Construct path for this specific run's version of the storm
                storm_run_path = path.join(MAXSS_working_directory, "output", run_name, 
                                           "maxss/storm-atlas/ibtracs", region, year_name, storm_name)
                
                fluxengine_files = sorted(glob(os.path.join(storm_run_path, "*.nc")))
                
                if not fluxengine_files:
                    print(f"   -> Run {run_name} not found. Skipping.")
                    continue

                # Print update to the console
                #print(f'Calculating Dynamic Flux for: {storm_name} ({year_name}) {run_name}...')

                hourly_fluxes_list = []
                time_list = []
                storm_total_tg = 0.0
            
                # Loop through each daily flux file
                for file_path in fluxengine_files:
                    with xr.open_dataset(file_path) as ds:
                        if grid_areas is None:
                            grid_areas = calculate_ellipsoidal_areas(ds, res=0.25)
                        
                        #  Loop through the 24 hours inside each file
                        for t_idx in range(len(ds.time)):
                            # Select the specific hour (isel = index selection)
                            ds_hour = ds.isel(time=t_idx)
                            
                            # Use .sel() to pick the land mask values that match the 
                            # moving window's current lat/lon coordinates.
                            current_water_map = water_proportion.sel(
                                latitude=ds_hour.latitude, 
                                longitude=ds_hour.longitude, 
                                method="nearest")
                            
                            # Calculate mass for THIS hour
                            # (g/m2/day) * area * (1/24) * ocean proportion in cell = grams in this hour
                            hourly_mass_map = ds_hour.OF * grid_areas * (1.0 / 24.0) * current_water_map
                            
                            # Sum across space for this specific hour
                            hour_total_grams = float(hourly_mass_map.sum(skipna=True))
                            
                            if not np.isnan(hour_total_grams):
                                hourly_tg = hour_total_grams / 1e12
                                storm_total_tg += hourly_tg # Accumulate total Tg
                                
                                # Append to lists for the NetCDF output
                                hourly_fluxes_list.append(hourly_tg)
                                
                                # Extract the exact timestamp for this hour
                                t_val = ds_hour.time.values
                                t_sec = (t_val - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
                                time_list.append(t_sec)

                #print(f"   -> Total Carbon: {storm_total_tg:.6f} Tg ({len(time_list)} timesteps)")
                
                # --- SAVE RUN-SPECIFIC NETCDF ---
                out_folder = path.join(netcdf_output_root, region, year_name)
                if not os.path.exists(out_folder): os.makedirs(out_folder)
                
                # Format: [STORM_ID]_[STORM_NAME]_[RUN_NAME].nc
                processedFilePath = os.path.join(out_folder, f"{storm_name}_{run_name}.nc")
    
                with Dataset(processedFilePath, 'w') as ncout:
                    ncout.createDimension("time", len(time_list))
                    var_time = ncout.createVariable("time", "i8", ("time",))
                    var_time.units = "seconds since 1970-01-01 00:00:00"
                    var_time[:] = np.array(time_list)
                    
                    var_hourly = ncout.createVariable("Hourly_flux", "f8", ("time",))
                    var_hourly.units = "Tg C hr-1"
                    var_hourly[:] = np.array(hourly_fluxes_list)
                    
                    var_total = ncout.createVariable("Total_flux", "f8")
                    var_total.units = "Tg C"
                    var_total.assignValue(storm_total_tg)
    
                # --- UPDATE MASTER SUMMARY FOR CSV ---
                storm_id = (storm_name, year_name, region)
                if storm_id not in master_summary:
                    master_summary[storm_id] = {}
                master_summary[storm_id][run_name] = storm_total_tg
    
                print(f"      - {run_name}: {storm_total_tg:.6f} Tg")

# Convert the dictionary into a list of dictionaries for Pandas
final_rows = []
for (storm, year, region), run_data in master_summary.items():
    # Start the row with the identifying metadata
    row = {
        'Storm': storm,
        'Year': year,
        'Region': region
    }
    # Add each run's total as a new column
    row.update(run_data)
    final_rows.append(row)

# Create DataFrame and save
if final_rows:
    df_final = pd.DataFrame(final_rows)
    
    # Optional: Ensure columns are in a specific order (Metadata then Runs)
    cols = ['Storm', 'Year', 'Region'] + runs
    # Only keep columns that actually ended up in the data
    cols = [c for c in cols if c in df_final.columns]
    df_final = df_final[cols]
    
    output_path = path.join(output_base, "storm_component_flux_summary.csv")
    
    # Ensure directory exists before saving
    if not os.path.exists(output_base):
        os.makedirs(output_base)
        
    df_final.to_csv(output_path, index=False)
    print(f"\n{'='*40}")
    print(f"SUCCESS: Summary saved to {output_path}")
    print(f"{'='*40}")
else:
    print("\nNo data was processed. Check your directory paths.")

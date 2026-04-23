# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 17:08:29 2026

@author: dw557
"""

import netCDF4 as nc
import numpy as np

def calculate_carbon_hourly_snapshot(file_path):
    """
    Calculates the mass of Carbon (Tg) for the first hour of the day.
    Logic: (Rate per day * Area) / 24 hours / 1e12
    """
    ds = nc.Dataset(file_path)
    
    # 1. Extract Variables for the first time slice [0]
    flux_00 = ds.variables['OF'][0, :, :] 
    lats = ds.variables['latitude'][:]
    lons = ds.variables['longitude'][:]
    
    # 2. Grid Area Calculation (m^2)
    R = 6371000  # Earth Radius
    dlat = np.abs(np.deg2rad(lats[1] - lats[0]))
    dlon = np.abs(np.deg2rad(lons[1] - lons[0]))
    
    lat_mesh, _ = np.meshgrid(lats, lons, indexing='ij')
    areas = (R**2) * np.cos(np.deg2rad(lat_mesh)) * dlat * dlon
    
    # 3. Calculate Mass for this specific hour
    # flux (g/m2/day) * area (m2) = grams/day rate
    total_grams_day_rate = np.nansum(flux_00 * areas)
    
    # Scale from "day rate" to "hourly mass"
    total_grams_hour = total_grams_day_rate / 24
    
    # Convert to Teragrams
    total_tg_hour = total_grams_hour / 1e12
    
    ds.close()
    return total_tg_hour

# --- File Paths ---
old_file = "E:/MAXSS_working_directory - 02_03_25 back up/output/MAXSS_RUN/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/MAXSS_2010176N16278_AL012010_ALEX_DOY_181_DATE_2010_06_30.nc"
new_file = "E:/MAXSS_working_directory/output/MAXSS_RUN/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/MAXSS_2010176N16278_AL012010_ALEX_DOY_181_DATE_2010_06_30.nc"

# --- Run ---
old_hourly = calculate_carbon_hourly_snapshot(old_file)
new_hourly = calculate_carbon_hourly_snapshot(new_file)

print(f"Carbon Mass (00:00-01:00) - Old File: {old_hourly:.8f} Tg")
print(f"Carbon Mass (00:00-01:00) - New File: {new_hourly:.8f} Tg")

# For your reference: The Total Daily Rate (what you had before)
print(f"\nDaily Rate (for comparison) - Old File: {old_hourly * 24:.8f} Tg/day")
print(f"Daily Rate (for comparison) - New File: {new_hourly * 24:.8f} Tg/day")
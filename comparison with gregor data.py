# -*- coding: utf-8 -*-
"""
Created on Fri Feb  6 10:28:12 2026

@author: dw557
"""

# Script to compare the results from my fluxengine run for Hurricane MARIA, with
# output from the gregor paper OceanSODAETHZv2

#Import required packages
import os
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime

# --- Paths ---
target_lon, target_lat = -73, 30
fluxengine_dir = "E:/MAXSS_working_directory/output/MAXSS_RUN/maxss/storm-atlas/ibtracs/north-atlantic/2017/2017260N12310_AL152017_MARIA"
oceansoda_flux_path = "E:/MAXSS_working_directory/OceanSODAETHZv2/fgco2-2010s-8D_25km-OceanSODAETHZv2.2024r01.nc"
oceansoda_dfco2_path = r"E:/MAXSS_working_directory/OceanSODAETHZv2/dfco2-2010s-8D_25km-OceanSODAETHZv2.2024r01.nc"
oceansoda_kw_path = r"E:/MAXSS_working_directory/OceanSODAETHZv2/kw-2010s-8D_25km-OceanSODAETHZv2.2024r01.nc"

# --- 1. Load FluxEngine Data ---
fluxengine_files = sorted([f for f in os.listdir(fluxengine_dir) if f.endswith('.nc')])

with nc.Dataset(os.path.join(fluxengine_dir, fluxengine_files[0]), 'r') as ds:
    fe_lons, fe_lats = ds['longitude'][:], ds['latitude'][:]
    lat_idx = (np.abs(fe_lats - target_lat)).argmin()
    lon_idx = (np.abs(fe_lons - target_lon)).argmin()

fe_flux, fe_dpco2, fe_kw, fe_time = [], [], [], []

for file in fluxengine_files:
    with nc.Dataset(os.path.join(fluxengine_dir, file), 'r') as ds:
        fe_flux.extend(ds.variables['OF'][:, lat_idx, lon_idx])
        fe_dpco2.extend(ds.variables['dpCO2'][:, lat_idx, lon_idx])
        fe_kw.extend(ds.variables['OK3'][:, lat_idx, lon_idx])
        
        u, c = ds.variables['time'].units, getattr(ds.variables['time'], 'calendar', 'standard')
        fe_time.extend(nc.num2date(ds.variables['time'][:], units=u, calendar=c))

fe_dates = [datetime.datetime(d.year, d.month, d.day, d.hour, d.minute) for d in fe_time]

# --- 2. Load OceanSODA Data (Generic loader) ---
def get_oceansoda_ts(path, var_name):
    with nc.Dataset(path, 'r') as ds:
        os_lats, os_lons = ds['lat'][:], ds['lon'][:]
        os_lat_i = (np.abs(os_lats - target_lat)).argmin()
        os_lon_i = (np.abs(os_lons - target_lon)).argmin()
        data = ds.variables[var_name][:, os_lat_i, os_lon_i]
        u, c = ds.variables['time'].units, getattr(ds.variables['time'], 'calendar', 'standard')
        dates = nc.num2date(ds.variables['time'][:], units=u, calendar=c)
        dates_fixed = [datetime.datetime(d.year, d.month, d.day, d.hour, d.minute) for d in dates]
        return dates_fixed, data

os_dates, os_flux_raw = get_oceansoda_ts(oceansoda_flux_path, 'fgco2')
_, os_dfco2 = get_oceansoda_ts(oceansoda_dfco2_path, 'dfco2')
_, os_kw = get_oceansoda_ts(oceansoda_kw_path, 'kw')

os_flux_gC = os_flux_raw * 0.012011 

# --- 3. Plotting ---
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 14), sharex=True)

#get lon/lat of comparison point for title
lat_str = f"{abs(target_lat)}°{'N' if target_lat >= 0 else 'S'}"
lon_str = f"{abs(target_lon)}°{'E' if target_lon >= 0 else 'W'}"

# Plot 1: Flux
# Use fe_dates (Python datetimes) and os_dates (Python datetimes)
ax1.plot(fe_dates, fe_flux, label='FluxEngine CO2 Flux', color='navy', lw=1.5)
ax1.step(os_dates, os_flux_gC, label='OceanSODA CO2 Flux', color='crimson', where='mid')

ax1.set_title(f'Hurricane Maria Comparison at {lat_str}, {lon_str}')

# Plot 2: Gradient
ax2.plot(fe_dates, fe_dpco2, label='FluxEngine dpCO2', color='teal', lw=1.5)
ax2.step(os_dates, os_dfco2, label='OceanSODA dfco2', color='darkorange', where='mid')

# Plot 3: Gas Transfer (kw)
ax3.plot(fe_dates, fe_kw, label='FluxEngine OK3 (N00)', color='darkgreen', lw=1.5)
ax3.step(os_dates, os_kw, label='OceanSODA kw (Fay et al.)', color='purple', where='mid')

# --- 4. Formatting the Axis ---
# Ensure we define the locator AFTER the plots are drawn
day_locator = mdates.DayLocator(interval=5)
major_formatter = mdates.DateFormatter('%Y-%m-%d')

ax3.xaxis.set_major_locator(day_locator)
ax3.xaxis.set_major_formatter(major_formatter)
ax3.xaxis.set_minor_locator(mdates.DayLocator(interval=1))

# Set the limits using datetime objects to match the axis type
ax3.set_xlim(fe_dates[0] - datetime.timedelta(days=1), 
             fe_dates[-1] + datetime.timedelta(days=1))

fig.autofmt_xdate(rotation=45)


for ax in [ax1, ax2, ax3]:
    # This brings the labels back for each panel
    ax.legend(loc='upper right', fontsize=10, frameon=True)
    
    # Ensure grid lines are visible on all plots
    ax.grid(True, which='major', linestyle='-', alpha=0.3)
    ax.grid(True, which='minor', linestyle=':', alpha=0.2)

plt.tight_layout()
plt.show()








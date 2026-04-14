# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 15:57:43 2026

@author: dw557
"""

## Script to plot total flux across storm

# Import required packages
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from pathlib import Path
import netCDF4 as nc

# Set working directory
MAXSS_working_directory = Path("E:/MAXSS_working_directory")
os.chdir(MAXSS_working_directory)

# Set save location for plots
output_dir_bar = Path("E:/MAXSS_working_directory/output/plots/total_flux_over_storm")
output_dir_bar.mkdir(parents=True, exist_ok=True)

output_dir_timeseries = Path("E:/MAXSS_working_directory/output/plots/flux_timeseries")
output_dir_timeseries.mkdir(parents=True, exist_ok=True)

# The list of runs to process
runs = ["MAXSS_RUN", "REF_RUN", "WIND_RUN", "SST_NO_GRADIENTS_RUN", 
        "SST_WITH_GRADIENTS_RUN", "SSS_RUN", "V_GAS_RUN", "PRESSURE_RUN"]

# Regions to create plots for
regions = ["north-atlantic"]

# Storms to skip
storms_to_skip = ["RINA", "MARIA"]

# Set up a colourblind friendly colour scheme
cb_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
run_colors = {
    "MAXSS_RUN": "k",                        # Black
    "REF_RUN": cb_colors[1],                 # Orange
    "WIND_RUN": cb_colors[7],                # Grey/Slate      
    "SST_NO_GRADIENTS_RUN": cb_colors[0],    # Blue 
    "SST_WITH_GRADIENTS_RUN": cb_colors[9],  # Light Blue/Cyan
    "SSS_RUN": cb_colors[4],                 # Purple
    "PRESSURE_RUN": cb_colors[5],            # Brown/Tan
    "V_GAS_RUN": cb_colors[2]}               # Green

# Set format for dates at bottom of plots
Month_Fmt = mdates.DateFormatter('%b %d')

def get_datetime(secondsSince1970):
    return datetime(1970, 1, 1) + timedelta(seconds=float(secondsSince1970))

# --- Reusable Plotting Function ---
def plot_flux_data(x_data, y_data, storm_name, ylabel, is_bar=True, output_path=None):
    """
    Handles both Bar charts (totals) and Line plots (timeseries).
    """
    fig, ax = plt.subplots(figsize=(8, 5) if is_bar else (8, 4))
    
    if is_bar:
        # Create Bar Chart
        bars = ax.bar(x_data, y_data, color='blue', width=0.8, edgecolor='white')
        for i, run_name in enumerate(x_data):
            if run_name in run_colors:
                bars[i].set_color(run_colors[run_name])
        
        ax.set_xticks(range(len(x_data)))
        ax.set_xticklabels(x_data, rotation=45, ha='right', fontsize=9)
        ax.set_xlabel('Model Run', fontsize=10)
    else:
        # Create Timeseries (y_data is a list of arrays)
        for i, run_name in enumerate(runs):
            clr = run_colors.get(run_name, 'k')
            ax.plot(x_data, y_data[i], color=clr, linewidth=1, label=run_name)
        
        # Zero reference line
        ax.plot(x_data, np.zeros(len(x_data)), "--", color="k", linewidth=1)
        ax.xaxis.set_major_formatter(Month_Fmt)
        ax.legend(fontsize=8, loc='best', ncol=2)
        ax.set_xlabel("Time", fontsize=10)

        # Set 2.5% custom padding 
        x_min, x_max = min(x_data), max(x_data)
        x_range = x_max - x_min
        padding = x_range * 0.025  # Calculate 2.5% of the total duration
        
        ax.set_xlim(x_min - padding, x_max + padding)        

        #set axis points
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=7))
        ax.xaxis.set_major_formatter(Month_Fmt) # Uses your '%b %d' format

    # Universal Styling
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(storm_name, fontsize=10)
    ax.tick_params(axis='both', which='major', labelsize=9, direction='in', length=2.5, width=1)
    ax.grid(True, linestyle=':', alpha=0.4)
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300)
    
    plt.show()
    plt.close(fig)

# --- Main Execution Block ---
if __name__ == "__main__":
    for region in regions:
        region_dir = MAXSS_working_directory / "maxss" / "storm-atlas" / "ibtracs" / region
        
        for year_folder in [d for d in region_dir.iterdir() if d.is_dir()]:
            year = year_folder.name
            for storm_folder in [s for s in year_folder.iterdir() if s.is_dir()]:
                storm = storm_folder.name
                
                if any(name in storm for name in storms_to_skip):
                    print(f"Skipping: {storm}")
                    continue
                
                print(f"Processing: {storm}")
                
                storm_totals = []
                storm_timeseries = []
                plot_dates = None

                for run_type in runs:
                    nc_file = Path(f"output/Spatially_integrated_fluxes/maxss/storm-atlas/ibtracs/{region}/{year}/{storm}_{run_type}.nc")
                    
                    if nc_file.exists():
                        with nc.Dataset(nc_file) as ds:
                            # Extract Data
                            total_val = float(ds.variables['Total_flux'][:])
                            hourly_flux = ds.variables['Hourly_flux'][:]
                            hourly_flux[hourly_flux == 0] = np.nan # Clean zeros
                            
                            storm_totals.append(total_val)
                            storm_timeseries.append(hourly_flux)
                            
                            if plot_dates is None:
                                plot_dates = [get_datetime(t) for t in ds.variables['time'][:]]
                    else:
                        storm_totals.append(0.0)
                        storm_timeseries.append(np.array([]))
                 
                # Set save locations 
                bar_save_path = output_dir_bar / f"Total_flux_over_storm_{storm}.png"
                timeseries_save_path = output_dir_timeseries / f"Total_flux_over_storm_{storm}.png"
                
                # 1. Plot Comparison Bar Chart
                
                plot_flux_data(runs, storm_totals, storm, 
                               ylabel='Total flux (Tg C)', is_bar=True,
                               output_path=bar_save_path)

                # 2. Plot Multi-run Timeseries
                if plot_dates is not None:
                    plot_flux_data(plot_dates, storm_timeseries, storm, 
                                   ylabel="Total flux (Tg C hr$^{-1}$)", is_bar=False,
                                   output_path=timeseries_save_path)
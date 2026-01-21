# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 12:25:51 2026

@author: dw557
"""

import numpy as np
import netCDF4 as nc

#Load in a single MAXSS data file to create the mask.
fluxengine_1_deg_1_day_file = 'E:/MAXSS_working_directory/output/MAXSS_RUN/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/1_degree_resolution_output/1deg_MAXSS_2010176N16278_AL012010_ALEX_DOY_160_DATE_2010_06_09.nc'

fluxengine_ds = nc.Dataset(fluxengine_1_deg_1_day_file)

#Load global data with xarray for easier masking
global_monthly_ds = nc.Dataset('E:/MAXSS_working_directory/global_background_all_3_months.nc')

global_ds_sept_flux = global_monthly_ds['flux'][:,:,0]

fluxengine_06_09_flux = fluxengine_ds['OF'][12,:,:]

#######################

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# --- Load Data (Your existing code) ---
# Load in a single MAXSS data file to create the mask.
fluxengine_1_deg_1_day_file = 'E:/MAXSS_working_directory/output/MAXSS_RUN/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/1_degree_resolution_output/1deg_MAXSS_2010176N16278_AL012010_ALEX_DOY_160_DATE_2010_06_09.nc'
fluxengine_ds = nc.Dataset(fluxengine_1_deg_1_day_file)

# Load global data with xarray for easier masking (Note: using netCDF4 as per your snippet)
global_monthly_ds = nc.Dataset('E:/MAXSS_working_directory/global_background_all_3_months.nc')

# Extract variables
global_ds_sept_flux = global_monthly_ds['flux'][:,:,0]
fluxengine_06_09_flux = fluxengine_ds['OF'][12,:,:]

# --- Plotting Code ---

# 1. Flatten the arrays to 1D
data_global = np.array(global_ds_sept_flux).flatten()
data_fluxengine = np.array(fluxengine_06_09_flux).flatten()

# 2. Remove NaN values (common in masked NetCDF data)
data_global = data_global[~np.isnan(data_global)]
data_fluxengine = data_fluxengine[~np.isnan(data_fluxengine)]

# 3. Create the Plot
plt.figure(figsize=(10, 6))

# Plot first histogram (Global)
plt.hist(data_global, bins=50, alpha=0.5, label='Global Sept Flux (Dan F)', color='blue', edgecolor='none')

# Plot second histogram (FluxEngine)
plt.hist(data_fluxengine, bins=50, alpha=0.5, label='FluxEngine 06/09 Flux', color='orange', edgecolor='none')

# Add labels and legend
plt.xlabel('Flux (g C m-2 day-1)')
plt.ylabel('Frequency')
plt.title('Comparison of Global Flux vs FluxEngine Flux across study area')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.3)

# Show plot
plt.show()
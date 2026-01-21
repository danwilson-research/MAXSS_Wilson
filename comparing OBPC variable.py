# -*- coding: utf-8 -*-
"""
Created on Fri Jan 16 14:06:34 2026

@author: dw557
"""

#Code to check for differences between flux engine model runs
#In this case between my FluxEngine output and that by Dan Ford
#this script works with data that has been regridded onto 1deg grid
#code created with assistance from Gemini AI

#Import required packages
import numpy as np
import netCDF4 as nc
import xarray as xr
import matplotlib.pyplot as plt
import os


#Resampled data
resampled_datafile = "E:/MAXSS_working_directory/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/Resampled_for_fluxengine_Ford_et_al_pco2.nc"

resampled_data = nc.Dataset(resampled_datafile)

plt.hist(resampled_data['pCO2water_mean'][0,:,:])

#load in fluxengine output
fluxengine_output = "E:/MAXSS_working_directory/output/MAXSS_RUN/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/MAXSS_2010176N16278_AL012010_ALEX_DOY_166_DATE_2010_06_15.nc"

fluxengine_output = nc.Dataset(fluxengine_output)

plt.hist(fluxengine_output['OBPC'][0,:,:])



import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# 1. Load the Datasets
resampled_path = "E:/MAXSS_working_directory/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/Resampled_for_fluxengine_Ford_et_al_pco2.nc"
fluxengine_path = "E:/MAXSS_working_directory/output/MAXSS_RUN/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/MAXSS_2010176N16278_AL012010_ALEX_DOY_166_DATE_2010_06_15.nc"

resampled_data = nc.Dataset(resampled_path)
fluxengine_output = nc.Dataset(fluxengine_path)

# 2. Extract and Flatten Data
# We access the variables, take the first time step [0,:,:], and flatten to 1D array
# We also filter out masked/NaN values to avoid plotting errors
pco2_data = resampled_data['pCO2water_mean'][0,:,:].flatten()
obpc_data = fluxengine_output['OBPC'][0,:,:].flatten()

# Remove masked values if the data is a MaskedArray
if np.ma.is_masked(pco2_data):
    pco2_data = pco2_data.compressed()
if np.ma.is_masked(obpc_data):
    obpc_data = obpc_data.compressed()

# 3. Create the Plot
plt.figure(figsize=(10, 6))

# Plot Resampled Data (Blue)
plt.hist(pco2_data, bins=50, alpha=0.5, label='Resampled Dan F OBPC var (input data)', color='blue', density=True)

# Plot FluxEngine Output (Orange)
plt.hist(obpc_data, bins=50, alpha=0.5, label='FluxEngine Output OBPC', color='orange', density=True)

# Formatting
plt.title('Comparison of Resampled pCO2 vs FluxEngine OBPC')
plt.xlabel('Value')
plt.ylabel('Frequency (Density)')
plt.legend(loc='upper right')
plt.grid(True, alpha=0.3)

plt.show()

# Close datasets
resampled_data.close()
fluxengine_output.close()



# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 14:30:42 2026

@author: dw557
"""

#Script to extract the air-sea co2 flux for the hurricane area from Dan F
#Global carbon budget data (monthly resolution)

#Load required packages
import os
import netCDF4 as nc
import datetime

## 1: extract the months data that we need

#Location of fluxengine data
fluxengine_data_path = "E:/MAXSS_working_directory/output/MAXSS_RUN/maxss/storm-atlas/ibtracs/north-atlantic/2010/2010176N16278_AL012010_ALEX/"

#Load in the FluxEngine flux data
fluxengine_files = os.listdir(fluxengine_data_path)

months_to_extract = []

# Loop through FluxEngine data
for file in fluxengine_files:
    one_day_of_data = nc.Dataset(fluxengine_data_path + file)

    # load in flux data
    #flux_data = one_day_of_data.variables['OF'][:]
    
    # load in one time from each day (in seconds since 1970,01,01)
    date = float(one_day_of_data.variables['time'][0])
    
    #get date from above float
    current_datetime = datetime.datetime.fromtimestamp(date, datetime.timezone.utc)

    days_since = (current_datetime - datetime.datetime(1970, 1, 15, tzinfo=datetime.timezone.utc)).total_seconds() / 86400
    

    #add to running list
    months_to_extract.append(current_datetime.month)
    
#extract the unique months and sort in order
unique_months = list(set(months_to_extract))
unique_months.sort()


## 2: extract the months of data from Dan F data

#location of Dan F global monthly flux data
global_monthly_data_path = 'E:/MAXSS_working_directory/Fordetal_UExP-FNN-U_surface-carbonate-system_v2025-1.nc'

#load in .nc of all monthly global data
global_monthly_flux = nc.Dataset(global_monthly_data_path)

#load in time variable
time_var = global_monthly_flux.variables['time']

needed_year_months = set((d.year, d.month) for d in all_dates)
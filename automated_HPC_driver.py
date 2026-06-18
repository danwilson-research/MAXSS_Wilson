#!/usr/bin/env python3
import os
import sys
from os import path

# 1. Catch the year and storm name passed from the custom SLURM batch file
if len(sys.argv) < 3:
    print("Error: Missing arguments. Usage: python3 driver.py [YEAR] [STORM_NAME]")
    sys.argv = ["2010", "COLIN"] # Safe fallback defaults for local debugging

passed_year = sys.argv[1]
passed_storm = sys.argv[2]

# 2. Pipeline Module Imports
import MAXSS_resample as re
import MAXSS_run as ru
import calculate_hourly_flux_across_storm as c_flux

# 3. Directory Routing Configurations
MAXSS_working_directory = "/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_working_directory/"
downloadedRoot = "/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_working_directory/Ford_et_al_GBC_fco2/flux"
configfiletemplate = "/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_Wilson/MAXSS_configuration_file_template.conf"
output_base = '/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_working_directory/output/Spatially_integrated_fluxes'
netcdf_output_root = path.join(MAXSS_working_directory, "output/Spatially_integrated_fluxes/maxss/storm-atlas/ibtracs")

runs = ["MAXSS_RUN", "REF_RUN", "WIND_RUN", "SST_NO_GRADIENTS_RUN",
        "SST_WITH_GRADIENTS_RUN", "SSS_RUN", "V_GAS_RUN", "PRESSURE_RUN"]
MAXSS_regions = ["north-atlantic"]

# 4. Overriding the internal settings with the active target storm parameters
specified_years = [passed_year]
specified_storms = [passed_storm]

verbose = True
test_run = False # Change to True if you want a quick 5-day verification run

print(f"\n[PIPELINE START]: Processing Year {passed_year} | Storm {passed_storm}")

# Step A: Data Resampling with automatic bypass check

# Extract active region name dynamically from your target list
active_region = MAXSS_regions[0]

# We pick a specific file that proves the resampling finished successfully last time
sentinel_file = path.join(
    MAXSS_working_directory, 
    f"maxss/storm-atlas/tropical/ibtracs/{active_region}/{passed_year}/{passed_storm}/Resampled_for_fluxengine_Ford_et_al_pco2_pre_storm_reference.nc" 
    # ^ NOTE: Double-check your exact internal resampling path/filenames inside MAXSS_resample if needed!
)

if os.path.exists(sentinel_file):
    print(f"[PIPELINE INFO]: Resampled data already exists for {passed_storm}. Skipping Step A to save time.")
else:
    print(f"[PIPELINE INFO]: Resampled data missing. Executing Step A (Data Resampling)...")
    re.MAXSS_resample_main(MAXSS_working_directory, downloadedRoot, specified_storms, MAXSS_regions, specified_years)

# Step B: Main FluxEngine Executions
ru.MAXSS_flux_run(MAXSS_working_directory, configfiletemplate, verbose, specified_storms, test_run)

# Note: Keeping hourly calculations commented out until you finish upgrading its internal script
# c_flux.calc_hourly_flux(MAXSS_working_directory, output_base, netcdf_output_root, runs, MAXSS_regions, storms_to_skip)

print(f"[PIPELINE END]: Successfully finished all jobs for {passed_storm}.\n")

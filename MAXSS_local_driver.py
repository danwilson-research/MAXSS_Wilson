# -*- coding: utf-8 -*-
"""

Driver script for running the MAXSS analysis locally

@Original author: Daniel Ford
@Editing author: Daniel Wilson
"""

import MAXSS_resample as re
import MAXSS_run as ru
import calculate_hourly_flux_across_storm as c_flux
from os import path

MAXSS_working_directory = "E:/MAXSS_working_directory"
downloadedRoot = "E:/MAXSS_working_directory/Ford_et_al_GBC_fco2/flux"
configfiletemplate="E:/MAXSS_Wilson/MAXSS_configuration_file_template.conf"
output_base = 'E:/MAXSS_working_directory/output/Spatially_integrated_fluxes'
netcdf_output_root = path.join(MAXSS_working_directory, "output/Spatially_integrated_fluxes/maxss/storm-atlas/tropical/ibtracs")
runs = ["MAXSS_RUN", "REF_RUN", "WIND_RUN", "SST_NO_GRADIENTS_RUN",
        "SST_WITH_GRADIENTS_RUN", "SSS_RUN", "V_GAS_RUN", "PRESSURE_RUN"]
MAXSS_regions = ["north-atlantic"]
specified_years = ['2010']

verbose = False

#Specify which storms you would like to run # if no storms specified, all storms run
specified_storms = [] #["RINA","AL052010_","ALEX" "BONNIE", "MARIA", , "COLIN", "DANIELLE"]

#When set to True, only MAXSS_main run is computed and only first 5 days modelled.
test_run = True

#re.MAXSS_resample_main(MAXSS_working_directory,downloadedRoot, specified_storms, MAXSS_regions, specified_years)

ru.MAXSS_flux_run(MAXSS_working_directory,configfiletemplate,verbose,specified_storms,test_run)

#c_flux.calc_hourly_flux(MAXSS_working_directory,output_base,netcdf_output_root,runs,MAXSS_regions,storms_to_skip)


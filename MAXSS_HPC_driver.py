"""

Driver script for running the MAXSS analysis locally

@Author: Daniel Wilson, based on teomplate from Daniel Ford

"""

# Import required packages and functions
import MAXSS_resample as re
import MAXSS_run as ru
import calculate_hourly_flux_across_storm as c_flux
from os import path

MAXSS_working_directory = "/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_working_directory/"
downloadedRoot = "/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_working_directory/Ford_et_al_GBC_fco2/flux"
configfiletemplate="/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_Wilson/MAXSS_configuration_file_template.conf"
output_base = '/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_working_directory/output/Spatially_integrated_fluxes'
netcdf_output_root = path.join(MAXSS_working_directory, "output/Spatially_integrated_fluxes/maxss/storm-atlas/tropical/ibtracs")
runs = ["MAXSS_RUN", "REF_RUN", "WIND_RUN", "SST_NO_GRADIENTS_RUN",
        "SST_WITH_GRADIENTS_RUN", "SSS_RUN", "V_GAS_RUN", "PRESSURE_RUN"]
MAXSS_regions = ["north-atlantic"]
specified_years = ['2010']

verbose = False

#Specify which storms you would like to run # if no storms specified, all storms run
specified_storms = ["ALEX"] #["RINA", "BONNIE", "MARIA", "ALEX", "COLIN", "AL052010_"]

#When set to True, only MAXSS_main run is computed and only first 5 days modelled.
test_run = True

## TESTING SCRIPT
import importlib.util

# Load the file directly from its exact filesystem location
spec = importlib.util.spec_from_file_location("resample_test", r"E:\MAXSS_Wilson\resample_test.py")
re = importlib.util.module_from_spec(spec)
spec.loader.exec_module(re)

print("Import successful! Available functions/vars:", dir(re))

#### remove testing section onceno longer testing 

re.MAXSS_resample_main(MAXSS_working_directory,downloadedRoot, specified_storms, MAXSS_regions,specified_years)

#ru.MAXSS_flux_run(MAXSS_working_directory,configfiletemplate,verbose,specified_storms,test_run)

#c_flux.calc_hourly_flux(MAXSS_working_directory, output_base, netcdf_output_root,runs, MAXSS_regions, specified_storms)


"""

Driver script for running the MAXSS analysis locally

@Author: Daniel Wilson, based on teomplate from Daniel Ford

"""

# Import required packages and functions
#import MAXSS_resample as re
import MAXSS_run_Dieke as ru
import calculate_hourly_flux_across_storm as c_flux
from os import path
import importlib.util
import sys

## Require explicit storm argument — exit with an error if missing
if len(sys.argv) < 2:
    print("ERROR: No storm name provided! Usage: python MAXSS_HPC_driver_Deike.py <STORM_NAME>")
    sys.exit(1)

specified_storms = [sys.argv[1]]
print(f"Running driver for storm: {specified_storms}")

## Manually import the updated resampling script.
# Path to your resample script
resample_script_path = r"/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_Wilson/MAXSS_resample_sig_wv_ht.py"

# Load the module directly from disk
spec = importlib.util.spec_from_file_location("re", resample_script_path)
if spec is None or spec.loader is None:
    raise FileNotFoundError(
        f"Could not find or load script at: {resample_script_path}")

re = importlib.util.module_from_spec(spec)
spec.loader.exec_module(re)

print("Successfully loaded MAXSS_resample_sig_wv_ht!")

MAXSS_working_directory = "/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_working_directory/"
downloadedRoot = "/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_working_directory/Ford_et_al_GBC_fco2/flux"
configfiletemplate="/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_Wilson/MAXSS_configuration_file_template_DEIKE.conf"
output_base = '/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_working_directory/output/Spatially_integrated_fluxes'
netcdf_output_root = path.join(MAXSS_working_directory, "output/Spatially_integrated_fluxes/maxss/storm-atlas/tropical/ibtracs")
runs = ["MAXSS_RUN", "REF_RUN", "WIND_RUN", "SST_NO_GRADIENTS_RUN",
        "SST_WITH_GRADIENTS_RUN", "SSS_RUN", "V_GAS_RUN", "PRESSURE_RUN"]
MAXSS_regions = ["north-atlantic"]
specified_years = ['2019']

swh_directory = "/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_working_directory/ERA5_sig_wv_ht/"

verbose = False

#Specify which storms you would like to run # if no storms specified, all storms run
#specified_storms = ["DORIAN"] #, "LORENZO", "HUMBERTO", "JERRY", "BARRY", "ANDREA", "AL032019"]

#When set to True, only MAXSS_main run is computed and only first 5 days modelled.
test_run = False

### TESTING SCRIPT ##
#
#import importlib.util
#from pathlib import Path
#
## Provide the actual file path to the module
#file_path = "/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_Wilson/resample_test.py"
#
## Pass BOTH module name AND file path
#spec = importlib.util.spec_from_file_location("resample_test", file_path)
#
#if spec is None or spec.loader is None:
#    raise FileNotFoundError(f"Failed to load spec from path: {file_path}. Double-check the file exists!")
#
#re = importlib.util.module_from_spec(spec)
#spec.loader.exec_module(re)
#
#print("Import successful! Available functions/vars:", dir(re))
#
##### remove testing section onceno longer testing ####

#re.MAXSS_resample_main(MAXSS_working_directory=MAXSS_working_directory,
#    downloadedRoot=downloadedRoot, swh_directory=swh_directory,  # <-- Pass swh_directory keyword argument
#    specified_storms=specified_storms, MAXSS_regions=MAXSS_regions, specified_years=specified_years,)

ru.MAXSS_flux_run(MAXSS_working_directory, configfiletemplate, verbose, specified_storms, test_run, rate_parameterisation='/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_Wilson/Deike_2018_parameterisation.py')

#c_flux.calc_hourly_flux(MAXSS_working_directory, output_base, netcdf_output_root,runs, MAXSS_regions, specified_storms)

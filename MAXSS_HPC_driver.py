import MAXSS_resample as re
import MAXSS_run as ru
import calculate_hourly_flux_across_storm as c_flux
from os import path

MAXSS_working_directory = "/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_working_directory/"
downloadedRoot = "/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_working_directory/Ford_et_al_GBC_fco2/flux"
configfiletemplate="/nobackup/beegfs/workspace/dw557/MAXSS_Wilson/MAXSS_configuration_file_template.conf"
output_base = '/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_working_directory/output/Spatially_integrated_fluxes'
netcdf_output_root = path.join(MAXSS_working_directory, "output/Spatially_integrated_fluxes/maxss/storm-atlas/ibtracs")
runs = ["MAXSS_RUN", "REF_RUN", "WIND_RUN", "SST_NO_GRADIENTS_RUN",
        "SST_WITH_GRADIENTS_RUN", "SSS_RUN", "V_GAS_RUN", "PRESSURE_RUN"]
MAXSS_regions = ["north-atlantic"]

verbose = True

#Specify which storms you would like to run # if no storms specified, all storms run
specified_storms = ["ALEX"] #["RINA", "BONNIE", "MARIA", "ALEX", "COLIN", "AL052010_"]

#When set to True, only MAXSS_main run is computed and only first 5 days modelled.
test_run = True

re.MAXSS_resample_main(MAXSS_working_directory,downloadedRoot, specified_storms, MAXSS_regions)

#ru.MAXSS_flux_run(MAXSS_working_directory,configfiletemplate,verbose,specified_storms,test_run)

#BUG TO FIX# hourly flux calculation script needs updating to include test run and use specified storms before it is used

#c_flux.calc_hourly_flux(MAXSS_working_directory,output_base,netcdf_output_root,runs,MAXSS_regions,storms_to_skip)


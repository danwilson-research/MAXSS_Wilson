#!/usr/bin/env python3
# This python script is to set up a series of SLURM SBATCH submissions automatically

import os
import time
import subprocess

## Step one - read the directory containing input data for fluxengine ##
MAXSS_regions = ["north-atlantic"]
specified_years = ["2011"]

MAXSS_working_directory = "/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_working_directory/"
SLURM_standard_script = "/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_Wilson/SLURM_standard_script.sh"

# Define the new central root directory for all slurm batch files
SLURM_FILES_ROOT = os.path.join(MAXSS_working_directory, "slurm_files")
SLURM_LOGS_ROOT = os.path.join(MAXSS_working_directory, "slurm_job_logs")

# Loop through the regions in MAXSS storm dataset
for region in MAXSS_regions:
    # Define the directory for the region
    region_directory = os.path.join(MAXSS_working_directory, "maxss", "storm-atlas", "tropical", "ibtracs", region)

    # Get a list of the years safely (and sorted)
    year_list = sorted([
        d for d in os.listdir(region_directory) 
        if os.path.isdir(os.path.join(region_directory, d))
    ])

    #### Loop through the years in the MAXSS storm dataset
    for year in year_list:
        
        # Filter based on specified years
        if len(specified_years) > 0 and year not in specified_years:
            print(f"Skipping year: {year}")
            continue
        
        current_year_dir = os.path.join(region_directory, year)
        
        # Get a list of the storms for THIS specific year (and sorted)
        storm_list = sorted([
            d for d in os.listdir(current_year_dir) 
            if os.path.isdir(os.path.join(current_year_dir, d))
        ])
        
        print(f"\n--- Region: {region} | Year: {year} ---")
        print("Storms found:", storm_list)

        # Loop through the storms to create custom SLURM scripts ----
        for storm in storm_list:
            
            # 1. Read the contents of your baseline SLURM template
            with open(SLURM_standard_script, 'r') as file:
                slurm_content = file.read()
            
            # 2. Target your exact lines and replace them with the dynamic storm name
            old_name_line = "#SBATCH --job-name=single-core_throttle_job          # A name for the job to be used when Job Monitoring "
            new_name_line = f"#SBATCH --job-name={storm}          # A name for the job to be used when Job Monitoring "
            
            # Define the absolute directory where you want the Slurm logs to go
            log_destination_dir = os.path.join(
                SLURM_LOGS_ROOT, "storm-atlas", "tropical", "ibtracs", region, year)
                
            os.makedirs(log_destination_dir, exist_ok=True) # Automatically generate the target folder tree
            
            old_out_line = "#SBATCH --output=test-single-core_throttle_-%j.out    # Standard output and error log"
            new_out_line = f"#SBATCH --output={log_destination_dir}/{storm}-%j.out    # Standard output and error log"
            
            old_echo_line = 'echo "Running a program on $SLURM_JOB_NODELIST"'
            new_echo_line = f'echo "Running {storm} on $SLURM_JOB_NODELIST"'
            
            slurm_content = slurm_content.replace(old_name_line, new_name_line)
            slurm_content = slurm_content.replace(old_out_line, new_out_line)
            slurm_content = slurm_content.replace(old_echo_line, new_echo_line)
            
            # 3. Inject the command line arguments for driver.py
            slurm_content = slurm_content.replace("__YEAR__", year)
            slurm_content = slurm_content.replace("__STORM__", storm)
            
            # 4. Define the dedicated slurm files directory for this region and year
            slurm_full_category_dir = os.path.join(
                SLURM_FILES_ROOT, "storm-atlas", "tropical", "ibtracs", region, year)
            
            os.makedirs(slurm_full_category_dir, exist_ok=True) 
            
            custom_script_path = os.path.join(slurm_full_category_dir, f"submit_{storm}.sh")
            
            # 6. Write the modified content out to the new file
            with open(custom_script_path, 'w') as file:
                file.write(slurm_content)
                
            print(f"   --> Generated custom SLURM script: {custom_script_path}")
            
            #7. AUTOMATIC SUBMISSION TO ATHENA CLUSTER QUEUE ...
            try:
                # 1. Execute 'sbatch' and capture both stdout and stderr
                result = subprocess.run(
                    ["sbatch", custom_script_path], 
                    capture_output=True, 
                    text=True, 
                    check=True
                )
                stdout_msg = result.stdout.strip()
                stderr_msg = result.stderr.strip()
                print(f"    [SUBMITTED]: {stdout_msg}")
                print(f"    [SUBMITTED]: {stderr_msg}")

#            except subprocess.CalledProcessError as e:
#                # If sbatch fails (exit code status non-zero), check=True raises this exception
#                print(f"    [SUBMISSION FAILED] Status code: {e.returncode}")
#                print(f"    [STDOUT]: {e.stdout.strip()}")
#                print(f"    [STDERR]: {e.stderr.strip()}")
                
                # 2. Extract the JobID using regex (matches the trailing digits)
                import re
                job_match = re.search(r"\d+$", stdout_msg)
                
                if job_match:
                    job_id = job_match.group()
                    
                    # Small grace pause: give the Slurm controller daemon a moment to register the state
                    time.sleep(2.5)
                    
                    # 3. Check 'squeue' explicitly for this JobID to confirm it's alive (Running or Pending)
                    verify_check = subprocess.run(
                        ["squeue", "-j", job_id], 
                        capture_output=True, 
                        text=True
                    )
                    
                    if job_id in verify_check.stdout:
                        print(f"    [VERIFIED]: Job {job_id} is safely tracking in the active queue.")
                    else:
                        # The job was submitted but isn't in squeue (it failed instantly)
                        print(f"    ? [CRITICAL ERROR]: Job {job_id} dropped out of the queue immediately after submission!")
                        
                        # Grab the accounting log to show you why it dropped
                        sacct_check = subprocess.run(
                            ["sacct", "-j", job_id, "--format=State,ExitCode", "-n", "-P"],
                            capture_output=True, text=True
                        )
                        print(f"    [DEAD STATE DETAILS]: {sacct_check.stdout.strip()}")
                else:
                    print("    ?? [WARNING]: Could not parse a valid JobID from Slurm string output.")
                    
            except subprocess.CalledProcessError as e:
                print(f"    ? [QUEUE FAILED]: Could not submit {storm}. Error: {e.stderr.strip()}")

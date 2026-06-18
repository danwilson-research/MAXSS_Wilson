#!/bin/bash -l
# Use bash and pickup a basic login environment.

#SBATCH --partition=hmq                    # Specifies that the job will run on the default queue nodes.
#SBATCH --job-name=single-core_throttle_job          # A name for the job to be used when Job Monitoring 
#SBATCH --time=48:00:00                     # maximum run time for a job hrs:min:sec
#SBATCH --nodes=1                           # Number of full nodes to use
#SBATCH --ntasks-per-node=1                 # Run a single task on each node
#SBATCH --ntasks=1                          # Run a single task (increase this value for parallelisation across CPUs)
#SBATCH --mem=23000                         # Memory per node specification is in MB. It is optional. 
#SBATCH --output=test-single-core_throttle_-%j.out    # Standard output and error log
#SBATCH --mail-type=ALL                     # Mail events (NONE)
#SBATCH --requeue                          # Specifies that the job will be requeued after a node failure.
                                            # The default is that the job will not be requeued.

pwd; hostname; date

echo "Running a program on $SLURM_JOB_NODELIST"
 
module load Workspace/v1

# Setup Conda
source ~/miniforge3/etc/profile.d/conda.sh
conda activate fluxengine_4_2

# Change to your target save location
cd /nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_Wilson/
 
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
 
#python test.py
python3 -u automated_HPC_driver.py __YEAR__ __STORM__

date
#!/bin/bash
#--------------------------------------------------------------------------------
#  SBATCH CONFIG
#--------------------------------------------------------------------------------
#SBATCH --job-name=1_convert_dicom    # name for the job
#SBATCH --cpus-per-task=4             # number of cores
#SBATCH --mem=64G                      # total memory
#SBATCH --time 0-00:15                # time limit in the form days-hours:minutes
#SBATCH --mail-user=mizzzousrl@gmail.com  # email address for notifications
#SBATCH --mail-type=END,FAIL

#SBATCH --partition General           # max of 1 node and 2 hours; use `Lewis` for larger jobs
#--------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"

## Module Commands

#path='/storage/hpc/group/sleeplab/raw'

## Run the matlab script
SCRIPT='1_convert_dicom.sh'

## Test the first 2 of pilot

#srun bash ${SCRIPT} SP021 Visit_1
#srun bash ${SCRIPT} SP022 Visit_1

#srun bash ${SCRIPT} SPIN2_013 V1

srun bash ${SCRIPT} SPIN2_016 V1
srun bash ${SCRIPT} SPIN2_023 V1

echo "### Ending at: $(date) ###


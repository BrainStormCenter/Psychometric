#!/bin/bash
#--------------------------------------------------------------------------------
#  SBATCH CONFIG
#--------------------------------------------------------------------------------
#SBATCH --job-name=2_preprocessing    # name for the job
#SBATCH --cpus-per-task=1             # number of cores
#SBATCH --mem=64G                     # total memory
#SBATCH --time 0-02:00                # time limit in the form days-hours:minutes
#SBATCH --mail-user=mizzzousrl@gmail.com  # email address for notifications
#SBATCH --mail-type=END,FAIL

#SBATCH --partition General           # max of 1 node and 2 hours; use `Lewis` for larger jobs
#SBATCH --licenses=matlab:1
#--------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"

## Module Commands
module load matlab/matlab-R2018b
module list

## Run the matlab script
SCRIPT=/storage/hpc/group/sleeplab/Psychometric/automated_processing/s02_preprocessing.m
DIRNAME=/storage/hpc/group/sleeplab/preprocessed/SPIN2_013_V1_unknown/

cd /storage/hpc/group/sleeplab/Psychometric/automated_processing/
srun matlab -nodesktop -nosplash -nodisplay -r "run('${SCRIPT} ${DIRNAME}');exit"
mv $DIRNAME ${DIRNAME}.smooth

echo "### Ending at: $(date) ###


#!/bin/bash
#SBATCH --job-name=Optimization_pleated_filter
#SBATCH --output=optimization.%j.out # %j expands to slurm JobID
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cluster=stheno
#SBATCH --partition=dms-cpu
#SBATCH --mem-per-cpu 16G
#load the correct modules
module load matlab/2020b
/afs/cad.njit.edu/sw.common/matlab-2020b/bin/matlab -nodisplay nosplash -r "Optimization_pleated_filter"

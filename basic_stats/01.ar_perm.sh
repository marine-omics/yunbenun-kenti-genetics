#!/bin/bash --login
#SBATCH --account=pawsey1088
#SBATCH --job-name=snpstats
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --array=0-999%200


module load singularity/4.1.0-slurm

mkdir -p perms

singularity run $MYSOFTWARE/sifs/ldne.sif Rscript run_permute.R ${SLURM_ARRAY_TASK_ID}


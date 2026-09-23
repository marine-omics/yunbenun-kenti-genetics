#!/bin/bash --login
#SBATCH --account=pawsey1088
#SBATCH --job-name=ldne
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:30:00
#SBATCH --array=1-150%200

module load singularity/4.1.0-slurm

mkdir -p jackknife_ds

args=$(sed "${SLURM_ARRAY_TASK_ID}q;d" grid_jk_ds.tsv)

singularity run $MYSOFTWARE/sifs/ldne.sif Rscript run_ne_jk_ds.R $args


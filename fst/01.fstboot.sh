#!/bin/bash --login
#SBATCH --account=pawsey1169
#SBATCH --job-name=fst
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=40G
#SBATCH --time=2:00:00

module load singularity/4.1.0-slurm

singularity run $MYSOFTWARE/sifs/ldne.sif Rscript fst_boot.R


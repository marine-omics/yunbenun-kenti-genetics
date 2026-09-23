#!/bin/bash --login
#SBATCH --account=pawsey1088
#SBATCH --job-name=snpstats_all
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --array=0-999%200

# Reruns the all-individuals permutations of 01.ar_perm.sh with the allelic
# richness rarefaction depth held fixed and matched to the kin-removed run.
# Requires both hfstat.rds and hfstat_nokin.rds. Results go to perms_all/ so
# the original perms/ are left untouched.

module load singularity/4.1.0-slurm

mkdir -p perms_all

singularity run $MYSOFTWARE/sifs/ldne.sif Rscript run_permute_all.R ${SLURM_ARRAY_TASK_ID}

#!/bin/bash

#SBATCH --job-name=cseqtl_model
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=12-00:00:00
#SBATCH --array=5-13382%250
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org
#SBATCH --output=Rout/model_%j_%A.out
#SBATCH --error=Rerr/model_%j_%A.err

# Load R module
module load fhR/4.2.0-foss-2021b

ARRAY="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/eQTL/eqtl_array_all_chr_rmdup.csv"

# Run the R script with the gene name and chromosome number
Rscript cseqtl_modelling.R $ARRAY $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK

#!/bin/bash

#SBATCH --job-name=trecase_model
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=2
#SBATCH --array=3-13368%700
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org
#SBATCH --output=Rout/out_trecase_%j_%A.out
#SBATCH --error=Rerr/out_trecase_%j_%A.err

# Load R module
module load fhR/4.2.0-foss-2021b

ARRAY="/fh/working/hsu_l/Mari/cseqtl/scripts/trecase_array_all_chr_rmdup.csv"

# Run the R script with the gene name and chromosome number
Rscript trecase_output.R $ARRAY $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK

#!/bin/bash
#SBATCH --job-name=agg_genes
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --output=Rout/agg_%j.out
#SBATCH --error=Rerr/agg_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjohnso5@fredhutch.org

# Load necessary modules
module load fhR/4.2.0-foss-2021b

# Run the R script
Rscript agg_genes_new.R


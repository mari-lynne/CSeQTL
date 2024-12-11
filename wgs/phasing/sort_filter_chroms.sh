#!/bin/bash

#SBATCH --job-name=Sort_Index
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=12
#SBATCH --array=1-22
#SBATCH --output=Rout/sort_%A_%a.out  
#SBATCH --error=Rerr/sort_%A_%a.err    
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

# Load necessary modules
ml BCFtools

# Set vars and directories
CHROM="chr${SLURM_ARRAY_TASK_ID}"
# THREADS=${SLURM_CPUS_PER_TASK}

# IN_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased/chunks"
OUT_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased"

TEMP_DIR="${OUT_DIR}/temp"
OUT_FILE="${OUT_DIR}/${CHROM}.phased.bcf" # .filtered

# Sort and index the filtered BCF
echo "Sorting and Indexing ${OUT_FILE}"
bcftools sort ${OUT_FILE} --temp-dir ${TEMP_DIR} --write-index
if [ $? -ne 0 ]; then
    echo "Error sorting ${OUT_FILE}"
    exit 1
fi

# Cleanup temporary files and dirs

echo "${CHROM} processing completed successfully."


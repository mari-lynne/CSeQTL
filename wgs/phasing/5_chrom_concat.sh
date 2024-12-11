#!/bin/bash

#SBATCH --job-name=concat_phased
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=12
#SBATCH --output=Rout/concat-%j-%a.out  
#SBATCH --error=Rerr/concat-%j-%a.err    
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

# Aims:
# Concatenate chromosomes

ml BCFtools

# 1) Set Directories and variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

THREADS=${SLURM_CPUS_PER_TASK}

IN_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased
OUT_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased
OUT_FILE=lls_sct_concat # _filtered
OUT_FILE2=lls_sct_concat_filtered


# Concatenate Chromsomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bcftools concat ${IN_DIR}/chr{1..22}.phased.bcf --threads ${THREADS} --naive -Ob -o ${OUT_DIR}/${OUT_FILE}.bcf
echo "Chromsome Concatenating - DONE!"
  if [ $? -ne 0 ]; then
      echo "Error processing chromosome ${i}"
      exit 1
  fi

bcftools index ${OUT_DIR}/${OUT_FILE}.bcf --threads ${THREADS}
echo "File Index - DONE!"
  if [ $? -ne 0 ]; then
      echo "Error processing Indexing BCF File"
      exit 1
  fi

# Also concatenate filtered file - used for PCA downstream

bcftools concat ${IN_DIR}/chr{1..22}.phased.filtered.bcf --threads ${THREADS} -Ob -o ${OUT_DIR}/${OUT_FILE2}.bcf
bcftools index ${OUT_DIR}/${OUT_FILE2}.bcf --threads ${THREADS}

echo "Chr Concatenating and Indexing - DONE!"
  
# Generate sample array file for next step
# bcftools query -l ${OUT_DIR}/${OUT_FILE}.bcf > ${OUT_DIR}/${OUT_FILE}.txt


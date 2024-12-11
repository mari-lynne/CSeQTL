#!/bin/bash

#SBATCH --job-name=concat_pca
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=20GB
#SBATCH --output=Rout/pca-%j-%A-%a.out  
#SBATCH --error=Rerr/pca-%j-%A-%a.err   
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

# Aims: Prep data PCA data then run

# Load required modules
module load PLINK2/20210826-linux_x86_64
module load BCFtools

# Variables and Directories
THREADS=${SLURM_CPUS_PER_TASK}

IN_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased
PCA_DIR=${IN_DIR}/pca
IN_FILE=lls_sct_concat_filtered

START_TIME=$(date)
echo "Beginning WGS PCA prep: ${START_TIME}"
echo ${IN_DIR}
echo ${PCA_DIR}
echo ${IN_DIR}/${IN_FILE}.bcf

# Ensure PCA_DIR exists
mkdir -p ${PCA_DIR}


# PCA Prep =====================================================================

# 1) Convert to bfile and remove duplicates
plink2 \
  --bcf ${IN_DIR}/${IN_FILE}.bcf \
  --rm-dup force-first \
  --set-all-var-ids @:#:\$r:\$a \
  --make-bed \
  --out ${IN_DIR}/${IN_FILE}_rmdup
  
# # Remove Duplicate variants
# plink2 \
#   --bfile ${IN_DIR}/${IN_FILE} \
#   --rm-dup force-first \
#   --make-bed \
#   --out ${IN_DIR}/${IN_FILE}_rmdup

# OPTION: Remove outliers after 1st pass and QC check in R
plink2 \
  --bfile ${IN_DIR}/${IN_FILE} \
  --remove /fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/metadata/geno_outliers.txt \
  --make-bed \
  --out ${IN_DIR}/${IN_FILE}_rmdup
  
# 2a) Calculate LD
plink2 \
  --bfile ${IN_DIR}/${IN_FILE}_rmdup \
  --indep-pairwise 50 5 0.2 \
  --out ${PCA_DIR}/${IN_FILE}_rmdup

# 2b) Prune variants  
plink2 \
  --bfile ${IN_DIR}/${IN_FILE}_rmdup \
  --extract ${PCA_DIR}/${IN_FILE}_rmdup.prune.in \
  --make-bed \
  --out ${PCA_DIR}/${IN_FILE}_rmdup_ld

END_TIME=$(date)
echo "Finished WGS PCA prep: ${END_TIME}"

# Run PCA ======================================================================

plink2 \
  --bfile ${PCA_DIR}/${IN_FILE}_rmdup_ld \
  --pca \
  --out ${PCA_DIR}/${IN_FILE}_pca_results

END_TIME=$(date)
echo "Finished WGS PCA: ${END_TIME}"

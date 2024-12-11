#!/bin/bash

# Aims:
# Filter participants from target study from the larger WHI cohort WGS data 
# This instance contains both SCT and LLS samples extracted from the now accessible TOPMed Freeze 12 WGS data (unphased)
# Also make a concatenated sample file for later uses in CSeQTL pipeline

#SBATCH --job-name=freeze12_sample_filter
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=20G
#SBATCH --output=Rout/%j-%a.out
#SBATCH --error=Rerr/%j-%a.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org
#SBATCH --array=1-23

# Load necessary modules
ml BCFtools

# 1) Set Directories and variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GENO_DIR=/fh/scratch/delete90/kooperberg_c/topmed_freeze12/minDP10 # BCF File Input Directory
GENO_FILE=freeze.12a.chr

IN_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results # Sample_file/manifest Directory
SAMPLE_FILE=${IN_DIR}/metadata/sct_lls_freeze12_ids.txt
OUT_DIR=${IN_DIR}/genotype/merged_study


# Check if directories exist
for DIR in "$GENO_DIR" "$IN_DIR" "$OUT_DIR"; do
  if [ ! -d "$DIR" ]; then
    echo "Directory $DIR does not exist. Exiting."
    exit 1
  fi
done

# 2) Filter participants across chromosomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
chr=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
chr_f=${chr[$SLURM_ARRAY_TASK_ID-1]} 
# Adjusted to zero-based indexing for bash, i.e chr1-1=0, the first element of our bash array

bcftools view -m2 -v snps --include 'FILTER="PASS" & INFO/AC>=2' --threads ${SLURM_CPUS_PER_TASK} --force-samples -S ${SAMPLE_FILE} ${GENO_DIR}/${GENO_FILE}${chr_f}.pass_and_fail.gtonly.minDP10.bcf --write-index -Ob -o ${OUT_DIR}/${GENO_FILE}${chr_f}_all_snps.bcf

# |\
# bcftools +fill-tags -Ob -- -t MAF | \
# bcftools view -i 'MAF>0.01' -Ob -o ${OUT_DIR}/${GENO_FILE}${chr_f}.bcf

if [ $? -ne 0 ]; then
  echo "Error in filtering chromosome ${chr_f}. Exiting."
  exit 1
fi

echo "Filtering for chromosome ${chr_f} - DONE!"

# rename chr x to chrom 23








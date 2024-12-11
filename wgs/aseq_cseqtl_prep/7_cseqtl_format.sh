#!/bin/bash
#SBATCH --job-name=cseqtl_snp_format
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=8
#SBATCH --time=4-0:00:00
#SBATCH --output=Rout/%j-%a.out  
#SBATCH --error=Rerr/%j-%a.err    
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org
#SBATCH --array=1-22

# Aims:
# Make recoded SNP genotype files per gene for CSeQTL

# Info:
# Array runs 22 instances of job, each corresponding to one chromosome
# Coordinate file is pre-generated in R (gene_coord_format.R) and subsetted by chromosome

# Set Directories and variables ================================================

ml BCFtools

IN_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased"
REF_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/resources"
OUT_DIR="${IN_DIR}/per_gene"
INDEX_FILE="${REF_DIR}/eqtl_500kb${SLURM_ARRAY_TASK_ID}.txt"

IN_FILE="${IN_DIR}/chr${SLURM_ARRAY_TASK_ID}.phased.filtered.bcf"
OUT_FILE="${IN_DIR}/chr${SLURM_ARRAY_TASK_ID}.phased.filtered.maf01.bcf"

# 1) Filter for MAF >= 0.01 ====================================================

# Apply MAF filter and write index
bcftools view ${IN_FILE} -q 0.01:minor -Oz --write-index -o ${OUT_FILE}

echo "MAF filtering - COMPLETE!"

# 3) Make CSeQTL SNP/Gene matrices =============================================

# Steps:
# For each line of the coordinate file, corresponding to a gene:
# Extract the chromosome, start, and end positions - use that to subset the filtered BCF file
# Then recode the genotype columns and save to a TSV file named after the gene ID

# NOTE: Check chromosome annotation: If BCF files use 'chr21', ensure to extract using 'chr${chr}'

while IFS=$'\t' read -r gene_id chr start end; do
    # Extract the variants for the given gene region
    bcftools view -r "chr${chr}:${start}-${end}" ${OUT_FILE} | \
    # Format the output directly into a tab-separated values file
    bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%GT\t%SAMPLE\n]' > ${OUT_DIR}/${gene_id}.tsv
done < ${INDEX_FILE}

echo "CSeQTL SNP per chromosome ${SLURM_ARRAY_TASK_ID} gene files - COMPLETE!"


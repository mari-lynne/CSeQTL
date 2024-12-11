#!/bin/bash
#SBATCH --job-name=snp_list_extraction   # Job name
#SBATCH --nodes=1                        # Number of nodes
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --cpus-per-task=4                # Number of CPU cores per task
#SBATCH --mem=8G                         # Memory per node
#SBATCH --time=2:00:00                   # Time limit hrs:min:sec
#SBATCH --output=log_%A_%a.out           # Standard output log
#SBATCH --error=log_%A_%a.err            # Standard error log
#SBATCH --array=1-22                     # Array range for chromosomes 1-22

# Load the BCFtools module
ml BCFtools

# Input directory
IN_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased"
# Output directory (adjust as necessary)
OUT_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased/snp_lists"

# Create the output directory if it doesn't exist
mkdir -p ${OUT_DIR}

# Specify input SNP file based on SLURM_ARRAY_TASK_ID
SNP_FILE="${IN_DIR}/chr${SLURM_ARRAY_TASK_ID}.phased.filtered.maf01.bcf"

# Specify output SNP list file
OUT_FILE="${OUT_DIR}/chr${SLURM_ARRAY_TASK_ID}_snp_list.txt"

# Extract SNP list using BCFtools
bcftools query -f '%CHROM:%POS:%REF:%ALT\n' ${SNP_FILE} > ${OUT_FILE}

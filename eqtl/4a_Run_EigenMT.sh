#!/bin/bash

#SBATCH --job-name=eigen_all
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-13382%250
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org
#SBATCH --output=log/eigen_all_cell%j_%A.out
#SBATCH --error=log/eigen_all_cell%j_%A.err

# Array lengths (see EigenMT_prep.R)

# Load required modules
ml fhR/4.2.0-foss-2021b
ml Python/3.8.2-GCCcore-9.3.0
ml scikit-learn/0.23.1-foss-2020a-Python-3.8.2

# Aims

# Per cell type run eigenMT to estimate the number of effective independent tests for downstream MT corrections

# Directories and Variables ----------------------------------------------------

array_file="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/eQTL/eqtl_array_all_chr_rmdup.csv" # eqtl_array_all_chr_rmdup.csv
# Extract gene_id and chromosome from the array file using the SLURM_ARRAY_TASK_ID, had to sub extra quotation marks
gene_id=$(awk -F, -v i="$SLURM_ARRAY_TASK_ID" 'NR == i+1 {gsub(/"/, "", $2); print $2}' "$array_file")
chromosome=$(awk -F, -v i="$SLURM_ARRAY_TASK_ID" 'NR == i+1 {gsub(/"/, "", $3); print $3}' "$array_file")

cell_type="All_snps"
# B_cells.csv  CD4_T_cells.csv  CD8_T_cells.csv  Monocytes_Macrophages.csv  Neutrophils.csv (from prep pt1)

# Input Data
snp_dir="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased/per_gene" # (formatted for CSeQTL input not yet eigenMT)

# Output Data
results_dir="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT"

echo "${SLURM_ARRAY_TASK_ID}"
echo "${gene_id}"
echo "${chromosome}"
echo "${cell_type}"

# # Run the R script -------------------------------------------------------------
# Rscript EigenMT_prep_pt2.R $gene_id $chromosome $cell_type $snp_dir $results_dir  # Add cell type variable

# Define file paths using the parsed variables
EIGEN_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/resources/eigenMT"
GEN_FILE="${results_dir}/temp/${gene_id}_genotype_matrix.txt"
GENPOS_FILE="${results_dir}/temp/${gene_id}_genotype_position.txt"
PHEPOS_FILE="${results_dir}/temp/${gene_id}_phenotype_position.txt"
QTL_FILE="${results_dir}/eigen_all_qtls.txt" # Update to include all results
OUT_FILE="${results_dir}/${cell_type}/${gene_id}_eigenMT_output.txt"
CHROM="chr${chromosome}"

# Would have to do a check in R in temp dir - if file exists read.csv if not then make it

#  # Run Eigen MT ----------------------------------------------------------------

# Check if all required files are present
if [[ -f "$GEN_FILE" && -f "$GENPOS_FILE" && -f "$PHEPOS_FILE" && -f "$QTL_FILE" ]]; then
    python "${EIGEN_DIR}/eigenMT.py" --CHROM $CHROM \
      --QTL $QTL_FILE \
      --GEN $GEN_FILE \
      --GENPOS $GENPOS_FILE \
      --PHEPOS $PHEPOS_FILE \
      --cis_dist 500000 \
      --window 50 \
      --OUT $OUT_FILE
else
    echo "Error: One or more required input files are missing."
    exit 1
fi

echo "Eigen MT complete for $cell_type SNPs"

# Clean files ------------------------------------------------------------------

# NA files are 31 bytes - Move all na files to  qtl_na folder

if [ -f "$OUT_FILE" ] && [ $(stat -c%s "$OUT_FILE") -eq 31 ]; then
    echo "Output file is 31 bytes, indicating NA results. Moving to qtl_na directory..."
    mkdir -p "${results_dir}/qtl_na"
    mv "$OUT_FILE" "${results_dir}/qtl_na/"
    echo "Moved $OUT_FILE to qtl_na."
fi

# Remove temp files

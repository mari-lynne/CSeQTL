#!/bin/bash

# Aims:
# - Fix header row: The first column (RS_Number) actually contains row numbers. We need to delete this column (excluding the header).
# - Add a column in the file to include the seed SNP ID extracted from the file name.
# - Filter for LD SNPs with a value greater than 0.75.

IN_DIR="/fh/working/hsu_l/Mari/ld_proxy/extra_snps"
OUT_DIR="${IN_DIR}/filtered"

# Create output directory if it doesn't exist
mkdir -p "${OUT_DIR}"

cd ${OUT_DIR}

# Loop through each text file ending with '_grch38.txt'
for file in ${IN_DIR}/*_grch38.txt; do
    # Extract the base name (Seed SNP rsID)
    base_name=$(basename "${file}" "_grch38.txt")

    # Process the file to adjust columns, add SEED_SNP, and filter by LD
    awk -v name="$base_name" 'BEGIN {FS=OFS="\t"} 
        NR == 1 {print $0, "SEED_SNP"}  # Print header as is, add new column name
        NR > 1 {  # For data rows
            print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, name  # Skip first column, add SEED_SNP, and print
        }' "$file" | awk '$7 >= 0.75' > "${OUT_DIR}/$(basename "${file%_grch38.txt}")_filtered.txt"
done

echo 'Files have been updated and filtered based on LD SNPs > 0.75.'

#!/bin/bash

# Aims:
# - Fix header row: The first column (RS_Number) actually contains row numbers. We need to delete this column (excluding the header).

IN_DIR="/fh/working/hsu_l/Mari/ld_proxy/filtered/AFR"
OUT_DIR="${IN_DIR}"

# Create output directory if it doesn't exist
mkdir -p "${OUT_DIR}"

cd ${OUT_DIR}

# Loop through all files ending with 'filtered.txt' in the specified directory
for file in ${IN_DIR}/*filtered.txt; do
# Process the file with AWK
awk 'BEGIN {FS=OFS="\t"}
        NR == 1 {  # Process the header
            $12 = "SEED_SNP"  # Rename the 12th column
            print  # Print the header
        }
        NR > 1 {  # Process the data rows
            # Start loop from the second column to skip the first column
            for (i = 2; i <= NF; i++) {
                printf "%s%s", $i, (i < NF ? OFS : ORS)  # Print fields with tab or newline
            }
        }
    ' "$file" > "${file%filtered.txt}2_filtered.txt"  # Output filename
done

echo 'All files have been processed.'

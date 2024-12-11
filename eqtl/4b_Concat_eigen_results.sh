#!/bin/bash

# Define the cell type and directories
cell_type="All_snps"
results_dir="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT"
input_dir="${results_dir}/${cell_type}"
output_file="${results_dir}/${cell_type}_eigen.txt"

# Get the header from the first file and write it to the output file
first_file=$(ls ${input_dir}/*.txt | head -n 1)
head -n 1 "$first_file" > "$output_file"

# Concatenate all files excluding their header
for file in "$input_dir"/*.txt; do
  tail -n +2 "$file" >> "$output_file"
done

echo "Consolidation complete. All files have been concatenated into $output_file"


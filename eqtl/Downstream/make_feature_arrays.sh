#!/bin/bash

# Set base directory
base_dir="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl"

# Define cell types
cell_types=("monocytes" "b_cells" "t_cells" "neuts")

# Loop through each cell type to generate array files
for cell_type in "${cell_types[@]}"; do

# Set the peak directory for the current cell type
peak_dir="${base_dir}/resources/chipseq/${cell_type}/lifted"

# List all .narrowPeak files and save them to the array file
array_file="${base_dir}/scripts/SCT/eqtl/Downstream/${cell_type}_array.txt"
ls "${peak_dir}"/*.narrowPeak > "$array_file"

num_features=$(wc -l < $array_file)
echo "Array file created for $cell_type: $array_file, Number of Features = ${num_features}"

done
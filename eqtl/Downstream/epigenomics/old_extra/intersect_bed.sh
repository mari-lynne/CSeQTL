#!/bin/bash

#SBATCH --job-name=intersect_bed
#SBATCH --time=1:00:00      
#SBATCH --mem=6G            
#SBATCH --cpus-per-task=1
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

# Aims:
# Extract SNPs overlapping with a CHIPseq feature

# Load BEDTools
ml BEDTools

# Set variables
cell_type="monocytes"

base_dir="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl"

peak_dir="${base_dir}/resources/chipseq/${cell_type}/lifted"

snp_dir="${base_dir}/results/eqtl/merged_study/snp_nexus/output/${cell_type}"

out_dir="${snp_dir}/peak_overlaps"

mkdir -p ${out_dir}
chmod -R u+w "${out_dir}"

# Get the full file extension of SNP bed file
bed_file=$(find "$snp_dir" -maxdepth 1 -type f -name "*.bed")

# Create an array of .narrowPeak BED files
peak_files=("${peak_dir}"/*.narrowPeak)

# Loop over all .narrowPeak files in peak_dir
for peak_file in "${peak_files[@]}"; do

  echo "Processing $bed_file with $peak_file"

  # Extract the basename of chipseq BED file (without extension)
  feature=$(basename "$peak_file" .narrowPeak)

  # Define the output file name
  out_file="${out_dir}/${cell_type}_${feature}_overlaps.txt"

  # Perform intersection and save output
  bedtools intersect -a ${bed_file} -b ${peak_file} -wa > ${out_file}
  echo "Intersection complete for $bed_file with $peak_file, results saved to ${out_file}."
  
done


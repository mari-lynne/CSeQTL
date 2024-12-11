#!/bin/bash

#SBATCH --job-name=intersect_permute
#SBATCH --time=1:00:00      
#SBATCH --mem=8G            
#SBATCH --cpus-per-task=1
#SBATCH --output=Rout/%x_%A_%a.out
#SBATCH --error=Rerr/%x_%A_%a.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org
#SBATCH --array=0-6  

# Load BEDTools
ml BEDTools

# Aims:
# Evaluate CSeQTL enrichment for CS-epigenetic features
# Count number of overlapping CS-eQTL SNPs with epi features
# Permute all eQTLs tested to get null distribution
# In R we will use x2 for significance testing

# run as sbatch enrich_permute2.sh cell_type

# Set variables ================================================================
cell_type=$1
base_dir="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl"
peak_dir="${base_dir}/resources/chipseq/${cell_type}/lifted"
snp_dir="${base_dir}/results/eqtl/merged_study/snp_nexus/output/${cell_type}"
out_dir="${snp_dir}/peak_overlaps"

# Create output directory if it doesn't exist
mkdir -p ${out_dir}
chmod -R u+w "${out_dir}" # grant write permissions

# Set file names ---------------------------------------------------------------

# INPUT:
# Get the epigenetic peak file using SLURM array
array_file="${base_dir}/scripts/SCT/eqtl/Downstream/${cell_type}_array.txt"
peak_file=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $array_file)
# Extract the feature name from the chipseq file
feature=$(basename "$peak_file" .narrowPeak)

# CSeQTL SNP BED file (Just run once)
eqtl_file=${snp_dir}/*.bed

# Background SNP BED file (Permuted)
snp_file="${base_dir}/results/eqtl/merged_study/snp_nexus/output/background_snps_38.bed"

# OUTPUT:
out_file="${out_dir}/eqtl_${cell_type}_${feature}_results.txt" # CSeQTL overlap results
# results_file="${out_dir}/${cell_type}_${feature}_results.csv" # overlap counts + permutation overlap counts

# Run CSeQTL SNP/Epi Intersections =============================================

echo "Processing $snp_file with $peak_file"

# Perform intersection and save output
bedtools intersect -a ${eqtl_file} -b ${peak_file} -wa -wb > ${out_file}
echo "Intersection complete for $snp_file with $peak_file, results saved to ${out_file}."

# Count the number of overlapping features and save to the overlap count results file
# overlap_count=$(wc -l < ${out_file})
# echo "$cell_type,${feature},0,${overlap_count}" >> $results_file

# TODO: Actually will need to redo to get n of overlapping features grouped by LD so we don't overcount

# # Permute Background SNP/Epi Intersections =====================================
# 
# # Number of iterations for permutation
# iterations=1000
# 
# # Get the number of SNPs in the original input so we are comparing same n of SNPs
# num_rows=$(wc -l < ${eqtl_file})
# 
# # Run 1000 iterations of permutation for the current feature
# for i in $(seq 1 $iterations); do
#   
#   temp_snp_file="${out_dir}/temp_snp_${i}.bed"
#   temp_out_file="${out_dir}/${cell_type}_${feature}_temp_intersect.txt"
#   
#   # Randomly subset the background SNP file with num_rows
#   shuf -n $num_rows $snp_file > $temp_snp_file
# 
#   # Run bedtools intersect with the subset SNPs
#   bedtools intersect -a ${temp_snp_file} -b ${peak_file} -wa > ${temp_out_file}
# 
#   # Count the number of overlapping SNPs
#   overlap_count=$(wc -l < ${temp_out_file})
# 
#   # Append the result to the results file
#   echo "$cell_type,${feature},${i},${overlap_count}" >> $results_file
# 
#   # Clean up temporary files
#   rm ${temp_snp_file} ${temp_out_file}
# done
# 
# echo "Processing complete for feature ${feature}"
# 
# # Print the final results file location
# echo "Results saved to $results_file"


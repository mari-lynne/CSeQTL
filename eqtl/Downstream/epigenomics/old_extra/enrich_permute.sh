#!/bin/bash

# Input files
snp_file="all_snps_19.bed"
epigenetic_file="E032-DNase.hotspot.broad.bed"

# Output results file
results_file="DNase_overlap_permute_results.txt"

# Number of rows to subset (n SNPs in orignal)
num_rows=2274

# Number of iterations
iterations=50

# Initialize the results file
echo "Iteration,OverlapCount" > $results_file

# Loop through the number of iterations
for i in $(seq 1 $iterations); do
  # Randomly subset the SNP file by num_rows and save to a temporary file
  temp_snp_file="all_snps_19_${i}.bed"
  shuf -n $num_rows $snp_file > $temp_snp_file

  # Run bedtools intersect and save the output to a temporary file
  bedtools intersect -a $temp_snp_file -b $epigenetic_file -wa -wb > DNase_overlap_permute_temp

  # Count the number of lines in the intersect file
  overlap_count=$(wc -l < DNase_overlap_permute_temp)

  # Append the result to the results file
  echo "${i},${overlap_count}" >> $results_file

  # Clean up the temporary SNP file and intersect file
  rm $temp_snp_file DNase_overlap_permute_temp
done

# Print the results file location
echo "Results saved to $results_file"

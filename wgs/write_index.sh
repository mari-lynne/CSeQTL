#!/bin/bash

#SBATCH --job-name=write_index
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-3033%200
#SBATCH --output=Rout/%j-%a.out  
#SBATCH --error=Rerr/%j-%a.err    
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

IN_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased/chunks"
THREADS=${SLURM_CPUS_PER_TASK}

ml BCFtools

INPUT_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${IN_DIR}/file_index.txt)
echo "Processing file: ${INPUT_FILE}"

# Print start time and date
echo "Job started at: $(date)"
start_time=$(date +%s)

# Write index ------------------------------------------------------------------

# Index bcf files
bcftools index ${IN_DIR}/${INPUT_FILE} --threads ${THREADS}


# Job Output -------------------------------------------------------------------

# Print end time and date
echo "Job ended at: $(date)"
end_time=$(date +%s)

# Calculate and print the duration
duration=$((end_time - start_time))
echo "Job duration: $duration seconds"

# Get and print memory usage using sacct
memory_usage=$(sacct --format=MaxRSS -j ${SLURM_JOB_ID} -n | awk '{print $1}')
echo "Memory usage for job ${SLURM_JOB_ID}: ${memory_usage}"

# Calculate memory efficiency
total_mem=$((${THREADS} * 20 * 1024)) # Total memory requested in MB
used_mem=$(echo $memory_usage | sed 's/K//' | awk '{print $1 / 1024}') # Convert KB to MB
efficiency=$(echo "scale=2; $used_mem / $total_mem * 100" | bc)
echo "Memory efficiency: $efficiency%"


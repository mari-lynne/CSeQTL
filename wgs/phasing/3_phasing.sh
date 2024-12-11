#!/bin/bash

#SBATCH --job-name=Phase_SCT_LLS
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=6
#SBATCH --array=1-3207%100 # Update to fit chunk_files.txt length # 3207%100
#SBATCH --output=Rout/phase_%j-%a.out  
#SBATCH --error=Rerr/phase_%j-%a.err    
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

# Note after 2_phasing_prep.sh run ls *.bcf > chunk_files.txt in chunks dir
# And cat chunk_files.txt | wc -l to get the number of chunk files which will be our max array value

# Setup vars and dirs ----------------------------------------------------------

THREADS=${SLURM_CPUS_PER_TASK}

EAGLE_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/resources/Eagle/Eagle_v2.4.1"
GENETIC_MAP="${EAGLE_DIR}/tables/genetic_map_hg38_withX.txt.gz"

IN_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/chunks"
REF_DIR="/fh/scratch/delete90/kooperberg_c/topmed_freeze10/phased"

OUT_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased"
mkdir -p ${OUT_DIR}


ml BCFtools
# Get the chunk file name from the file list and extract Chromosome
CHUNK_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${IN_DIR}/chunk_files.txt) # Extracts i'th line from file list
CHUNK_NAME=$(basename ${CHUNK_FILE} .bcf) # Extracts the chunk name without the bcf extension
CHR=$(echo ${CHUNK_NAME} | cut -d'_' -f1 | cut -d'r' -f2) # -d'_' -f1 = cuts the first field before the _delimiter, cut -d'r' -f2 using r as the delimter, extract the field after the delim 
OUT_NAME="${OUT_DIR}/chunks/${CHUNK_NAME}_phased"


# Extract start and end positions
CHUNK_START=$(echo ${CHUNK_NAME} | cut -d'_' -f2 | cut -d'-' -f1)
CHUNK_END=$(echo ${CHUNK_NAME} | cut -d'_' -f2 | cut -d'-' -f2)
# Define flanking size
FLANKING_SIZE=100000

echo "Starting Phasing of:"
echo "Chunk ${CHUNK_FILE}"
echo "Chunk Name ${CHUNK_NAME}"
echo "Chromosome ${CHR}"
echo "Output file = ${OUT_NAME}"

# Create output directories if they don't exist
mkdir -p ${OUT_DIR}
mkdir -p ${OUT_DIR}/chunks

cd ${EAGLE_DIR}

# Run Eagle to phase the chunk
./eagle --geneticMapFile ${GENETIC_MAP} \
      --vcfRef ${REF_DIR}/freeze.10b.chr${CHR}.pass_only.phased.bcf \
      --vcfTarget ${IN_DIR}/${CHUNK_FILE} \
      --allowRefAltSwap \
      --outPrefix ${OUT_NAME} \
      --numThreads ${THREADS} \
      --chrom chr${CHR} \
      --bpStart ${CHUNK_START} \
      --bpEnd ${CHUNK_END} \
      --bpFlanking ${FLANKING_SIZE} \
      --vcfOutFormat b 

# Index chunk
bcftools index ${OUT_NAME}.bcf --threads ${THREADS}


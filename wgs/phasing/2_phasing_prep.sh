#!/bin/bash

#SBATCH --job-name=phasing_prep
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-22
#SBATCH --output=Rout/%j-%a.out  
#SBATCH --error=Rerr/%j-%a.err    
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

ml BCFtools

# Define the input and output directories
IN_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study"
OUT_DIR="${IN_DIR}/chunks"
CHR=${SLURM_ARRAY_TASK_ID}
THREADS=${SLURM_CPUS_PER_TASK}

# Create an output directory for the chunks if it doesn't exist
mkdir -p ${OUT_DIR}

# Define the input and output file names
INPUT_FILE="${IN_DIR}/freeze.12a.chr${CHR}_all_snps.bcf"

# Define chunk size and overlap
CHUNK_SIZE=1000000  # 1MB
OVERLAP=100000  # 100kb

# Get chromosome length
CHR_LENGTH=$(bcftools view -h ${INPUT_FILE} | grep -P "^##contig=<ID=chr${CHR}," | sed -n 's/.*length=\([0-9]*\).*/\1/p')

echo "Chromosome ${CHR} Length = ${CHR_LENGTH} BP"

# Split the BCF file into chunks with overlap
START=1
while [ $START -lt $CHR_LENGTH ]; do
  END=$((START + CHUNK_SIZE - 1))
  if [ $END -gt $CHR_LENGTH ]; then
    END=$CHR_LENGTH
  fi
  REGION="chr${CHR}:${START}-${END}"
  CHUNK_OUTPUT="${OUT_DIR}/chr${CHR}_${START}-${END}.bcf"
  
  # Extract the region
  bcftools view -r ${REGION} ${INPUT_FILE} -Ob -o ${CHUNK_OUTPUT}

  # Update start position for next chunk with overlap
  START=$((START + CHUNK_SIZE - OVERLAP))
done

echo "BCF Splitting DONE"

# Note after 2_phasing_prep.sh run ls *.bcf > chunk_files.txt in chunks dir
# And cat chunk_files.txt | wc -l to get the number of chunk files which will be our max array value for phasing script

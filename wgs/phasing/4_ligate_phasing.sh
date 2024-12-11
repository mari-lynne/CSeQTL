#!/bin/bash

#SBATCH --job-name=Ligate_BCF_Chunks
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=6
#SBATCH --array=1-22
#SBATCH --output=Rout/ligate_%j.out  
#SBATCH --error=Rerr/ligate_%j.err    
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

# Load necessary modules
ml BCFtools

# Set directories
IN_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased/chunks" # Output from 3_phasing.sh
OUT_DIR="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased" # Where to save chunks concatenated into chromosomes (filtered/sorted/indexed)

TEMP_RESIZE_DIR="${IN_DIR}/resized_${CHROM}" # For resizing phased chunks
TEMP_SORT_DIR="${OUT_DIR}/temp" # For sorting process

# Make output dirs if they don't exist
mkdir -p ${OUT_DIR}
mkdir -p ${TEMP_RESIZE_DIR}
mkdir -p ${TEMP_SORT_DIR}

# Set variables and file names
CHROM="chr${SLURM_ARRAY_TASK_ID}" 
THREADS=${SLURM_CPUS_PER_TASK}
OVERLAP=100000 # Chunk phasing overlap size in bp (set in 3_phasing.sh)

OUT_FILE1="${OUT_DIR}/${CHROM}.phased"
OUT_FILE2="${OUT_DIR}/${CHROM}.phased.filtered"

# Adjust chunk windows ---------------------------------------------------------
adjusted_chunks=()
prev_start=0
prev_end=0

# Create a list of BCF chunks for the chromosome
CHUNK_FILES=($(ls ${IN_DIR}/${CHROM}*phased*.bcf | sort -V))

for file in "${CHUNK_FILES[@]}"; do
    # Extract start and end from filename
    filename=${file##*/}
    start=$(echo $filename | grep -oP '(?<=_)\d+(?=-)')
    end=$(echo $filename | grep -oP '(?<=-)\d+(?=_phased)')

    # Convert start and end to integers
    start=$((start))
    end=$((end))

    # Echo file name, start, end, and adjusted values
    echo "Processing file: $filename"
    echo "Start: $start"
    echo "End: $end"

    # Adjust start and end for overlaps
    if [ "$prev_start" -eq 0 ]; then
        new_start=$start
    else
        new_start=$((start + OVERLAP / 2))
    fi

    if [ "$prev_end" -ne 0 ]; then
        new_end=$((end - OVERLAP / 2))
    else
        new_end=$end
    fi

    echo "New Start: $new_start"
    echo "New End: $new_end"
    echo "Previous Start: $prev_start"
    echo "Previous End: $prev_end"

    adjusted_chunks+=("${CHROM}:${new_start}-${new_end}")
    prev_start=$start
    prev_end=$end
done

# Resize chunks ----------------------------------------------------------------

for i in "${!CHUNK_FILES[@]}"; do
    adjusted_region=${adjusted_chunks[$i]}
    output_file="${TEMP_RESIZE_DIR}/${CHUNK_FILES[$i]##*/}"
    echo "Resizing ${CHUNK_FILES[$i]} to ${output_file} with region ${adjusted_region}"
    bcftools view -r "$adjusted_region" -Ob -o "$output_file" "${CHUNK_FILES[$i]}"
    if [ $? -ne 0 ]; then
        echo "Error resizing chunk ${CHUNK_FILES[$i]}"
        exit 1
    fi
done

# Concatenate the resized chunks -----------------------------------------------
echo "Concatenating resized chunks for ${CHROM}"
RESIZED_FILES=($(ls ${TEMP_RESIZE_DIR}/${CHROM}_*_phased.bcf | sort -V))
bcftools concat "${RESIZED_FILES[@]}" --threads ${THREADS} -Ob -o ${OUT_FILE1}.bcf
if [ $? -ne 0 ]; then
    echo "Error concatenating chunks for ${CHROM}"
    exit 1
fi

# Filter the concatenated BCF --------------------------------------------------
echo "Filtering ${OUT_FILE1}"
bcftools view ${OUT_FILE1}.bcf --include 'FILTER="PASS" & INFO/AC>=2' --threads ${THREADS} -Ob -o ${OUT_FILE2}.bcf
if [ $? -ne 0 ]; then
    echo "Error filtering ${OUT_FILE1}"
    exit 1
fi

# Sort and index the filtered BCF
echo "Sorting and Indexing ${OUT_FILE2}"
bcftools sort ${OUT_FILE2}.bcf --temp-dir ${TEMP_SORT_DIR} -Ob -o ${OUT_FILE2}.bcf --write-index
if [ $? -ne 0 ]; then
    echo "Error sorting ${OUT_FILE2}"
    exit 1
fi

echo "Ligation of chunks for ${CHROM} completed. Output written to ${OUT_FILE}"

# Clean up temporary resized chunks
rm -rf ${TEMP_RESIZE_DIR}
rm -rf ${TEMP_SORT_DIR}





#!/bin/bash

# Set variables
cell_type="neuts"
refdir="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/resources"
outdir="${refdir}/chipseq/${cell_type}"

# In terminal, download liftover tools and chain files if needed
# make program executable chmod +x liftOver
# wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# run script in script dir where liftover program is saved ./liftover.sh

# Create an output directory for the converted files
lifted_outdir="${outdir}/lifted"
mkdir -p "${lifted_outdir}"  # Create directory if it doesn't exist

# 2) liftBed: LiftOver the bed file to a new build
for i in "${outdir}"/*.narrowPeak; do
  # Extract the base name of the file without the extension
  base_name=$(basename "${i}" .narrowPeak)

  # Define output file paths
  lifted_file="${lifted_outdir}/${base_name}.hg38.narrowPeak"
  unmapped_file="${lifted_outdir}/${base_name}.unmapped"

  # Run liftOver
  ./liftOver "${i}" "${refdir}/hg19ToHg38.over.chain.gz" -bedPlus=6 "${lifted_file}" "${unmapped_file}"

  # Check if the liftOver was successful
  if [ -s "${lifted_file}" ]; then
    echo "LiftOver completed successfully for ${i}."
  else
    echo "LiftOver failed for ${i}. Check ${unmapped_file} for details."
  fi
done

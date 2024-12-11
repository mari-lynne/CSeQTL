#!/bin/bash

# Set variables
refdir="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/resources"
in_dir="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/fuma"

# Create an output directory for the converted files
lifted_outdir="${in_dir}/lifted"
mkdir -p "${lifted_outdir}"  # Create directory if it doesn't exist

# In terminal, download liftover tools and chain files if needed
# make program executable chmod +x liftOver
# wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg38ToHg19.over.chain.gz

# run script in script dir where liftover program is saved ./liftover.sh

# Define output file paths
lifted_file="${lifted_outdir}/eQTL_snps_19.bed"
unmapped_file="${in_dir}/eQTL_snps_38.bed"

cd ${refdir} # has liftover program downloaded

# Run liftOver
./liftOver "${i}" "${refdir}/hg38ToHg19.over.chain.gz" "${lifted_file}" "${unmapped_file}"

# Check if the liftOver was successful
if [ -s "${lifted_file}" ]; then
echo "LiftOver completed successfully for ${i}."
else
  echo "LiftOver failed for ${i}. Check ${unmapped_file} for details."
fi
done

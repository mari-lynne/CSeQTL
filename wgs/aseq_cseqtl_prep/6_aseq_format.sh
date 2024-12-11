#!/bin/bash

#SBATCH --job-name=Aseq_Format
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=2
#SBATCH --output=Rout/Format-Aseq-%j-%a.out  
#SBATCH --error=Rerr/Format-Aseq-%j-%a.err 
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org
#SBATCH --array=1-2111%350   # Number of gwas samples

# Aims:
# Split bcf file into individual sample files
# Filter for heterozygous snps and reformat genotype column (|=tab delim)

# Info:
# Run in parallel jobs using slurm and sample array

# 1) Variables and Directories ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ml BCFtools # /1.14-GCC-11.2.0

IN_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased
IN_FILE=lls_sct_concat_filtered
SAMPLE_FILE=${IN_DIR}/lls_sct_samples.txt

OUT_DIR=${IN_DIR}/hap
mkdir -p ${OUT_DIR}

# 2) Assign samples from txt file  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

f=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_FILE})
echo "Sample ID = ${f}"

# Check if Sample ID is empty and exit if so
if [ -z "$f" ]; then
  echo "Error: Sample ID is empty. Exiting script."
  exit 1
fi

# 3) Make ASeq heterozygous SNP file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1) Subset the bcf file by sample $(f)
# 2) Query relevant genotype fields [TGT]
# 3) Filter phased heterozygous genotype information (-i'GT="0|1" | GT="1|0"')
# 4) Reformat genotype field to tab delim and save output files

# Start time
start_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Task started at: $start_time"

bcftools view -s ${f} ${IN_DIR}/${IN_FILE}.bcf \
      | bcftools query -f '%CHROM\t%POS[\t%TGT\n]\n' -i 'GT="0|1" | GT="1|0"' \
      | sed 's:|:\t:g' - > ${OUT_DIR}/${f}_hap.txt

echo "Hetrozygous SNP File (per sample) for ASeq - DONE!"

end_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Task ended at: $end_time"

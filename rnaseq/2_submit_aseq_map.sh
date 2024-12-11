#!/bin/bash

#SBATCH --array=1-809%250
#SBATCH --job-name=SCT_Map_ASeq
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12G
#SBATCH --output=Rout/sct-map-%j-%a.out
#SBATCH --error=Rerr/sct-map-%j-%a.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

ml fhR/4.2.2-foss-2021b

# Input
MANIFEST=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/metadata/sct_manifest.csv 
# Manifest contains 6 cols: index, geno_id, bam_id, rna_batch, geno_dir, bam_dir
REF_DIR=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/resources/ # Where exon_rds file is stored

# Output
OUT_DIR1=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/rnaseq/ # rnaseq mapped reads
OUT_DIR2=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/rnaseq/ASE/ # ASE mapped reads

mkdir -p ${OUT_DIR1}
mkdir -p ${OUT_DIR2}

echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Memory per CPU: $SLURM_MEM_PER_CPU"

Rscript 2_aseq_map.R ${MANIFEST} ${REF_DIR} ${OUT_DIR1} ${OUT_DIR2} ${SLURM_ARRAY_TASK_ID} ${SLURM_CPUS_PER_TASK} ${SLURM_MEM_PER_CPU}






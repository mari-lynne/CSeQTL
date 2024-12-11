#!/usr/bin/Rscript

### Packages ------------------------------------------------------------------
# Installed extra packages to rhino home dir https://sciwiki.fredhutch.org/rModules
library(stringr)
library(data.table)
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(CSeQTL)
library(asSeq)

###  Directories ---------------------------------------------------------------

# Inherit bash args
args <- commandArgs(trailingOnly = TRUE)
cat("Inherited arguments:", paste(args, collapse = " "), "\n")
cat("array ", Sys.getenv('SLURM_ARRAY_TASK_ID'), "\n")
cat("array ", args[5], "\n")
cat("job ", Sys.getenv('SLURM_JOB_NAME'), "\n")

# SLURM VARS
cpus <- as.integer(args[6])
mem_per_cpu <- as.integer(args[7])

# Read in Manifest and extract file names --------------------------------------
manifest <- read.csv(args[1], stringsAsFactors = FALSE)
thisRun <- manifest[as.integer(args[5]), ] # extract ith row of manifest csv

#### Input ---------------------------------------------------------------------

# IDs
geno_id <- as.character(thisRun$geno_id)
bam_id <- as.character(thisRun$bam_id)

# Input Directories
gen_dir <- as.character(thisRun$geno_dir)
bam_dir <- as.character(thisRun$bam_dir)
ref_dir <- args[2]

# Study
study <- as.character(thisRun$rna_batch)

# RNAseq Input File Name
if (study == "LLS") {
  bam_file <- paste0(bam_dir, bam_id, ".Aligned.sortedByCoord.out.md.bam")
} else {
  bam_file <- paste0(bam_dir, bam_id, ".accepted_hits.merged.markeddups.recal.bam")
}

# Exon Files # Made in 1_get-exon-info 
gtf_rds_fn <- paste0(ref_dir, "exon_by_genes.rds")
genes <- readRDS(gtf_rds_fn)

# Hap file
geno_file <- paste0(gen_dir, geno_id, "_hap.txt") # Made in 6_aseq_format.sh

#### Output --------------------------------------------------------------------

# Output Directories
rna_outdir <- args[3]
ase_outdir <- args[4]

# Output File Names
bam_filtered <- paste0(rna_outdir, bam_id, "-output.filtered.asSeq.bam")
bam_sort <- paste0(rna_outdir, bam_id, "-output.filtered.asSeq.sortQ")

# ASE Output
ase_out_file <- paste0(ase_outdir, geno_id, "-output.trecase.txt")

# Check names
cat("rna files: ", bam_file, bam_id, bam_filtered, bam_sort, "\n")
cat("wgs files: ", geno_id, geno_file, ase_out_file, "\n")

# 1) Filter RNA bam files for dups/low qual reads ------------------------------
start_time <- Sys.time()
start <- format(start_time, format = "%H:%M:%S")
print(paste("Starting ASE Pipeline:", start))

## Set up filter parameters
PE <- TRUE

flag1 <- Rsamtools::scanBamFlag(
  isUnmappedQuery = FALSE,
  isSecondaryAlignment = FALSE,
  isDuplicate = FALSE,
  isNotPassingQualityControls = FALSE,
  isSupplementaryAlignment = FALSE,
  isProperPair = PE
)

param1 <- Rsamtools::ScanBamParam(flag = flag1, what = "seq", mapqFilter = 255)

max_mem <- cpus * mem_per_cpu
max_mem <- max_mem * 0.8

## Filter bam file
Rsamtools::filterBam(bam_file,
                     destination = bam_filtered,
                     param = param1,
                     maxMemory = max_mem,
                     nThreads = cpus) # set in slurm script

end_time <- Sys.time()
time_taken <- format(end_time - start_time, format = "%H:%M:%S")
print(paste("bam filtering - DONE!", "Time taken:", time_taken))

# 2) Get total read count per gene ---------------------------------------------
start_time <- Sys.time()

bamfile <- Rsamtools::BamFileList(bam_filtered, yieldSize = 1000000)
se <- GenomicAlignments::summarizeOverlaps(
  features = genes,
  reads = bamfile,
  mode = "Union",
  singleEnd = !PE,
  ignore.strand = TRUE,
  fragments = PE
)

TReC <- as.data.frame(SummarizedExperiment::assay(se))
write.table(TReC, file = paste0(ase_outdir, geno_id, "-TReC.txt"), quote = FALSE, sep = "\t")

end_time <- Sys.time()
time_taken <- format(end_time - start_time, format = "%H:%M:%S")
print(paste("Total Read Count - DONE!", "Time taken:", time_taken))

# 3) Sort reads by Qname -------------------------------------------------------
start_time <- Sys.time()

sortBam(file = bam_filtered,
        destination = bam_sort,
        byQname = TRUE,
        maxMemory = max_mem,
        nThreads = cpus)

end_time <- Sys.time()
time_taken <- format(end_time - start_time, format = "%H:%M:%S")
print(paste("bam sorting - DONE!", "Time taken:", time_taken))

# 4) Get allele-specific read counts -------------------------------------------
start_time <- Sys.time()

asSeq::extractAsReads(
  input = paste0(bam_sort, ".bam"),
  snpList = geno_file,
  min.avgQ = 20,
  min.snpQ = 20
)

end_time <- Sys.time()
time_taken <- format(end_time - start_time, format = "%H:%M:%S")
print(paste("ASE counting - DONE!", "Time taken:", time_taken))

## 4b) Count allele-specific read counts ASRC ----------------------------------
start_time <- Sys.time()

se1 <- GenomicAlignments::summarizeOverlaps(
  features = genes,
  reads = sprintf("%s_hap1.bam", bam_sort),
  mode = "Union",
  singleEnd = !PE,
  ignore.strand = TRUE,
  fragments = PE
)
se2 <- GenomicAlignments::summarizeOverlaps(
  features = genes,
  reads = sprintf("%s_hap2.bam", bam_sort),
  mode = "Union",
  singleEnd = !PE,
  ignore.strand = TRUE,
  fragments = PE
)
seN <- GenomicAlignments::summarizeOverlaps(
  features = genes,
  reads = sprintf("%s_hapN.bam", bam_sort),
  mode = "Union",
  singleEnd = !PE,
  ignore.strand = TRUE,
  fragments = PE
)

hap1 <- as.data.frame(SummarizedExperiment::assay(se1))
hap2 <- as.data.frame(SummarizedExperiment::assay(se2))
hapN <- as.data.frame(SummarizedExperiment::assay(seN))
cts <- cbind(TReC, hap1, hap2, hapN)
dim(cts)
cts[1:2, ] # Check count file

write.table(cts, file = ase_out_file, quote = FALSE, sep = "\t", eol = "\n")

end_time <- Sys.time()
time_taken <- format(end_time - start_time, format = "%H:%M:%S")
print(paste("ASE mapping - DONE!", "Time taken:", time_taken))

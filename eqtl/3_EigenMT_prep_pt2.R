# Aims:

# Format data for input into Eigen MT
# Use Eigen MT to estimate number of effective tests
# Run per gene, as indicated by slurm array

library(data.table)
library(stringr)
library(dplyr)
library(tidylog)

## Read Args -------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
gene_id <- args[1]
chrom <- as.character(args[2])
cell_type <- args[3]
snp_dir <- args[4]
results_dir <- args[5]

# Remove any extra quotes from vars
print(gene_id)
print(chrom)

# chrom <- array[614, "chromosome_name"]
# gene_id <- array[614, "ensembl_gene_id"]

# Check if output files already exist
genotype_matrix_file <- file.path(results_dir, paste0("temp/", gene_id, "_genotype_matrix.txt"))
genotype_position_file <- file.path(results_dir, paste0("temp/", gene_id, "_genotype_position.txt"))
phenotype_position_file <- file.path(results_dir, paste0("temp/", gene_id, "_phenotype_position.txt"))

# If any of the files do not exist, proceed with the processing
if (!file.exists(genotype_matrix_file) || 
    !file.exists(genotype_position_file) || 
    !file.exists(phenotype_position_file)) {

  # Set number of threads for data.table
  setDTthreads(0L)
  
  # Read Genotype Data ----------------------------------------------------------
  GEN <- fread(file.path(snp_dir, paste0(gene_id, ".tsv")))
  colnames(GEN) <- c("chrom", "pos", "ref", "alt", "genotype", "sample_id")

  # Recode Genotypes to CSeQTL codes
  recode_genotype <- function(genotype) {
    ifelse(genotype == "0|0", 0,
           ifelse(genotype == "0|1", 1,
                  ifelse(genotype == "1|0", 1,
                         ifelse(genotype == "1|1", 2, NA))))
  }

  # Recode the genotype column
  set(GEN, j = "genotype", value = recode_genotype(GEN$genotype))
  set(GEN, j = "snp_id", value = paste(GEN$chrom, GEN$pos, GEN$ref, GEN$alt, sep = ":"))

  # Create GENPOS file -----------------------------------------------------------
  GENPOS <- select(GEN, snp_id, chrom, pos) %>% distinct()

  # Transpose GEN file for eigenMT input
  GEN <- pivot_wider(GEN, id_cols = snp_id, names_from = sample_id, values_from = genotype)

  # Read and process Phenotype Position File -------------------------------------
  
# Coords for CSeQTL were based on +/- 500MB; eigenMT needs original values
phe_file <- paste0("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/resources/eqtl_500kb", chrom, ".txt")
  PHEPOS <- fread(phe_file, header = FALSE)
  colnames(PHEPOS) <- c("gene_name", "chrom_probe", "s1", "s2")
  PHEPOS <- filter(PHEPOS, gene_name == gene_id)
  PHEPOS <- PHEPOS %>% mutate(s1 = (s1 + 500000),
                              s2 = s2 - 500000)
  PHEPOS$chrom_probe <- paste0("chr", PHEPOS$chrom_probe) 
  # Save files -------------------------------------------------------------------
  # Save Genotype Matrix
  fwrite(GEN, file = genotype_matrix_file, sep = "\t", quote = FALSE, row.names = FALSE)

  # Save Genotype (variant) Position File
  fwrite(GENPOS, file = genotype_position_file, sep = "\t", quote = FALSE, row.names = FALSE)

  # Save Phenotype (gene) position file
  fwrite(PHEPOS, file = phenotype_position_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
} else {
  print(paste0("Files for gene ", gene_id, " already exist. Skipping genotype and phenotype processing."))
}

# Filter eQTL results file ----------------------------------------------------
QTL <- read.csv(file = file.path(results_dir, paste0(cell_type, ".csv"))) %>% dplyr::rename(`p-value` = p.value)
QTL <- filter(QTL, gene_name == gene_id)

# Save Filtered QTL Results
fwrite(QTL, file = file.path(results_dir, paste0("temp/", cell_type, "_", gene_id, "_qtl_results.txt")),
       sep = "\t", quote = FALSE, row.names = FALSE)

# Test eigen -------------------------------------------------------------------

# python eigenMT.py --CHROM chr11 \
# --QTL "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT/temp/ENSG00000166311_qtl_results.txt" \
# --GEN "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT//temp/ENSG00000166311_genotype_matrix.txt" \
# --GENPOS "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT/temp/ENSG00000166311_genotype_position.txt" \
# --PHEPOS "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT/temp/ENSG00000166311_phenotype_position.txt" \
# --cis_dist 500000 \
# --window 50 \
# --OUT "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT/ENSG00000166311_eigenMT_output.txt"


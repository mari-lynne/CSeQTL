
library("ggplot2")
library("ggpubr")
library("viridis")
library("stringr")
library("data.table")
library("dplyr")
library("tidylog")

# ASE Results  -----------------------------------------------------------------

base_dir <- "/fh/working/hsu_l/Mari/cseqtl"
in_dir <- file.path(base_dir, "eqtl/trecase_only")
snp_dir <- file.path(base_dir, "genotype/phased/per_gene")
meta_dir <- file.path(base_dir, "metadata") 
setwd(in_dir)

# Inherit arguments from submission script
args <- commandArgs(trailingOnly = TRUE)
cat("Inherited arguments:", paste(args, collapse = " "), "\n")

array_file <- args[1]
index <- as.integer(args[2])
array <- read.csv(file = array_file)
gene_id <- array[index, "ensembl_gene_id"]
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))

# # For test
# array <- read.csv(file = "/fh/working/hsu_l/Mari/cseqtl/scripts/trecase_array_all_chr_rmdup.csv")
# index <- 1
# gene_id <- array[index, "ensembl_gene_id"]

### Load ASE data --------------------------------------------------------------

print(gene_id)
ase_out <- fread(sprintf("%s_eqtl.txt", gene_id))
# Marker ID doesnt have SNP ID in it :( checked input for trecase and rownames were snp ids
# Get original file - remake input

### Reload SNP data ------------------------------------------------------------

setDTthreads(0L)
SNP <- fread(paste0(snp_dir,"/", gene_id, ".tsv"))
colnames(SNP) <- c("chrom", "pos", "ref", "alt", "genotype", "sample_id")

# Recode bases to CSeQTL codes
recode_genotype <- function(genotype) {
  ifelse(genotype == "0|0", 0,
         ifelse(genotype == "0|1", 1,
                ifelse(genotype == "1|0", 3,
                       ifelse(genotype == "1|1", 4, 5))))
} # NOTE: Different codes used compared with CSeQTL package

SNP[, `:=`(genotype_code = recode_genotype(genotype),
           snp_id = paste(chrom, pos, ref, alt, sep = ":"))]

# Select the required columns, filter rows to include correct samples, then remove duplicate snps
ids <- read.csv(file = file.path(meta_dir, "cseqtl_ids_rna_order.csv"))
SNP <- SNP[sample_id %in% ids$geno_id, .(sample_id, snp_id, genotype_code)]
SNP <- SNP[, .(genotype_code = genotype_code[1]), by = .(snp_id, sample_id)]

# Reshape data to wide format using dcast
SNP <- dcast(SNP, snp_id ~ sample_id, value.var = "genotype_code")

# Convert to matrix for CSeQTL and update row names
r_names <- SNP[[1]]
SNP <- t(apply(SNP[,-1], 1, as.numeric))
rownames(SNP) <- r_names

### Add SNP IDs to ASE results -------------------------------------------------

# Filter for marker ID rows
SNP_2 <- SNP[ase_out$MarkerRowID, ]

# Add SNP ID to ase
ase_out$snp_id <- row.names(SNP_2)
ase_out <- ase_out %>% select(-GeneRowID)

# Rewrite new results file
write.table(ase_out, file = sprintf("%s_eqtl.txt", gene_id), sep = "\t", quote = F, row.names = F)

# Need to do this for all results, form a final results file
# Run 1_cseqtl_model_gene_submit.sh for edditing all files

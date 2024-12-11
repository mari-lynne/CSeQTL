# Aims:
# Read in genotype and aseq data
# Final formatting for CseQTL model
# Run CSeQTL
# Compute and save summary statistics

library("stringr")
library("data.table")
library("dplyr")
library("smarter")
library("CSeQTL")

# Inherit arguments
args <- commandArgs(trailingOnly = TRUE)
cat("Inherited arguments:", paste(args, collapse = " "), "\n")

array_file <- args[1]
index <- as.integer(args[2])
array <- read.csv(file = array_file)
gene_id <- array[index, "ensembl_gene_id"]
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))

# Data Directories -------------------------------------------------------------

base_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results"

snp_dir <- file.path(base_dir, "genotype/merged_study/phased/per_gene")

aseq_dir <- file.path(base_dir, "/genotype/merged_study/rnaseq/ASE/qc_genes/per_gene")

meta_dir <- file.path(base_dir, "metadata/merged")

out_dir <- file.path(base_dir, "/genotype/merged_study/eQTL")

# Load Covariates --------------------------------------------------------------

# See covar_modelling.R
ids <- read.csv(file = file.path(meta_dir, "cseqtl_ids_rna_order.csv"))
# model_covars_base
covars <- read.csv(file = file.path(meta_dir, "cseqtl_sample_xx_vars.csv"))
covars <- rename(covars, log_lib_size = lib.size)
num_vars <- covars %>% select(draw_age, log_lib_size)

# center numeric covariates (poss dont need to recentre pcs tho)
# Function to center a numeric variable
centre <- function(x) {
  if(is.numeric(x)) {
    return(x - mean(x, na.rm = TRUE))
  } else {
    return(x)
  }
}

num_vars <- as.data.frame(sapply(num_vars, centre))
covars <- covars %>% select(-draw_age, -log_lib_size)
model_covars <- cbind(num_vars, covars)
row.names(model_covars) <- ids$geno_id

# Load RHO (Cell type proportions) ---------------------------------------------
RHO <- read.csv(file.path(meta_dir, "cseqtl_rho.csv"))
row.names(RHO) <- ids$geno_id

# Combine T cells
# RHO <- RHO %>%
#   mutate(T_cells = CD4_T_cells + CD8_T_cells) %>%
#   select(-CD4_T_cells, -CD8_T_cells)

RHO <- as.matrix(RHO)

# Load SNP data ----------------------------------------------------------------

# Genotype data made in aggregate, gathers SNPs per gene
setDTthreads(0L)
SNP <- fread(paste0(snp_dir,"/", gene_id, ".tsv"))
colnames(SNP) <- c("chrom", "pos", "ref", "alt", "genotype", "sample_id")

recode_genotype <- function(genotype) {
  ifelse(genotype == "0|0", 0,
         ifelse(genotype == "0|1", 1,
                ifelse(genotype == "1|0", 2,
                       ifelse(genotype == "1|1", 3, 5))))
}

# Recode bases to CSeQTL codes
SNP[, `:=`(
  genotype_code = recode_genotype(genotype),
  snp_id = paste(chrom, pos, ref, alt, sep = ":")
)]

# Select the required columns and filter rows
SNP <- SNP[sample_id %in% ids$geno_id, .(sample_id, snp_id, genotype_code)]

### Filter duplicate SNPs ------------------------------------------------------
# Checked test gene (5000 SNPs 178 dup sites - majority were the same genotype i.e 0:0 or 3:3 a few 1/2s)
# Based on that it's fine to just keep the first entry (shouldn't majorly impact data)
SNP <- SNP[, .(genotype_code = genotype_code[1]), by = .(snp_id, sample_id)]

# Reshape data to wide format using dcast
SNP <- dcast(SNP, snp_id ~ sample_id, value.var = "genotype_code")

# Convert to matrix for CSeQTL and update row names
SNP <- as.matrix(SNP)
r_names <- SNP[,1]
SNP <- apply(SNP, 2, as.numeric)
row.names(SNP) <- r_names
SNP <- SNP[,-1]

# Load ASeq data ---------------------------------------------------------------

# Aseq data made in Aseq format > aggregate_genes.R
ase_file <- list.files(aseq_dir, pattern = gene_id)
aseq_dat <- fread(file.path(aseq_dir, ase_file))
colnames(aseq_dat) <- c("sample_id", "trec", "hap1", "hap2", "hapN")
aseq_dat$ASREC <- aseq_dat$hap1 + aseq_dat$hap2

## Check orders ----------------------------------------------------------------

# Order aseq data, covars, and RHO matrix by SNP data
aseq_dat <- aseq_dat[match(colnames(SNP), aseq_dat$sample_id), ]
model_covars <- model_covars[match(colnames(SNP), row.names(model_covars)), ]
RHO <- RHO[match(colnames(SNP), row.names(RHO)), ]

## Design Matrix --------------------------------------------------------------

# convert factor variables
model_covars$sct <- as.factor(model_covars$sct)
model_covars$ethnicity <- as.factor(model_covars$ethnicity)
model_covars$rna_batch <- as.factor(model_covars$rna_batch)
design <- model.matrix(~ ., data = model_covars)

# Make extra aseq vectors for model --------------------------------------------

TREC <-  aseq_dat$trec
names(TREC) <- aseq_dat$sample_id

ASREC <- aseq_dat$ASREC
names(ASREC) <- aseq_dat$sample_id

PHASE <- rep(1, nrow(aseq_dat))
names(PHASE) <- aseq_dat$sample_id

hap2 <- aseq_dat$hap2
names(hap2) <- aseq_dat$sample_id

# Run eQTL Model ---------------------------------------------------------------

cat("Starting CSeQTL Modelling for Gene:", gene_id)
start_time <- Sys.time()

# Model
gs_out <- CSeQTL_GS(XX=design, SNP = SNP, PHASE = PHASE,
                    TREC = TREC, hap2 = hap2, ASREC = ASREC, RHO = RHO,
                    ncores = n_cores, show = TRUE)

end_time <- Sys.time()
time_taken <- format(end_time - start_time, format = "%H:%M:%S")
print(paste("CSeQTL Modelling - COMPLETE!","Time taken:", time_taken))

# Extract Results --------------------------------------------------------------


# MJ additions Jan-2024
cistrans        = 0.01 # p-value cutoff to determine if the eQTLs is cis or trans
celltypes       = colnames(RHO)
snp_ids         = rownames(SNP)

MU_A <- gs_out$MU_A # Read count Allele A
MU_B <- gs_out$MU_B # Read count Allele B
ETA <- ifelse(MU_A == MU_B, 1, MU_B / MU_A) # Relative expression of Allele B/ Allele A

# Calculate P-values using chisq comparing model assumptions
PVAL_trec 			= pchisq(gs_out$LRT_trec[,,drop=FALSE], df = 1, lower.tail = TRUE)
PVAL_trecase 		= pchisq(gs_out$LRT_trecase[,,drop=FALSE], df = 1, lower.tail = TRUE)
PVAL_cistrans 	= pchisq(gs_out$LRT_cistrans[,,drop=FALSE], df = 1, lower.tail = TRUE)
CIS_eqtl 				= ifelse(PVAL_cistrans > cistrans, "CIS","TRANS") # 0.01
LRT_eqtl				= ifelse(CIS_eqtl == "CIS", gs_out$LRT_trecase, gs_out$LRT_trec) 
PVAL_eqtl				= pchisq(LRT_eqtl, df = 1, lower.tail = TRUE)

# Create a data-frame with significant results
results_mix <- data.frame(
  SNP = snp_ids,
  eQTL = CIS_eqtl, 
  MUa = MU_A,
  MUb = MU_B,
  ETA = ETA,
  LRT_eqtl = LRT_eqtl,
  PVAL_final = PVAL_eqtl
) # Using if else version based on cis/trans

results_ase <- data.frame(
  SNP = snp_ids,
  eQTL = CIS_eqtl, 
  MUa = MU_A,
  MUb = MU_B,
  ETA = ETA,
  LRT_ase = gs_out$LRT_trecase,
  PVAL_ase = PVAL_trecase
)

results_trec <-  data.frame(
  SNP = snp_ids,
  eQTL = CIS_eqtl, 
  LRT_trec = gs_out$LRT_trec,
  PVAL_trec = PVAL_trec
)


### Filter sig eQTLs -----------------------------------------------------------

# Filter SNPs with p-values below 0.05
sig_results <- results_mix[results_mix$PVAL_final.Neutrophils < 0.05 |
                             results_mix$PVAL_final.Monocytes_Macrophages < 0.05 |
                             results_mix$PVAL_final.CD4_T_cells < 0.05 |
                             results_mix$PVAL_final.CD8_T_cells < 0.05 |
                             results_mix$PVAL_final.B_cells < 0.05, ]

# TODO - check trecase only results as didn't save
sig_ase <- results_ase[results_ase$PVAL_ase.Neutrophils < 0.05 |
                             results_ase$PVAL_final.Monocytes_Macrophages < 0.05 |
                             results_ase$PVAL_final.CD4_T_cells < 0.05 |
                             results_ase$PVAL_final.CD8_T_cells < 0.05 |
                             results_ase$PVAL_final.B_cells < 0.05, ]

sig_trec <- results_trec[results_trec$PVAL_trec.Neutrophils < 0.05 |
                         results_trec$PVAL_trec.Monocytes_Macrophages < 0.05 |
                         results_trec$PVAL_trec.CD4_T_cells < 0.05 |
                         results_trec$PVAL_trec.CD8_T_cells < 0.05 |
                         results_trec$PVAL_trec.B_cells < 0.05, ]

write.table(sig_results, file = file.path(out_dir, paste0("combined/", gene_id, "_pvals.txt")),
            quote = F, sep = "\t", row.names = F)

write.table(sig_ase, file = file.path(out_dir, paste0("trecase/", gene_id, "_pvals.txt")),
            quote = F, sep = "\t", row.names = F)

write.table(sig_trec, file = file.path(out_dir, paste0("trec/", gene_id, "_pvals.txt")),
            quote = F, sep = "\t", row.names = F)

# Save Beta values
beta <- gs_out$PARS[row.names(gs_out$PARS) %in% sig_results$SNP, ]

write.table(beta, file = file.path(out_dir, paste0("combined/", gene_id, "_betas.txt")),
            quote = F, sep = "\t", row.names = F)


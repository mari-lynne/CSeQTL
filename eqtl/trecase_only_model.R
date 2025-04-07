# Linear model testing ---------------------------------------------------------

library(asSeq)
?asSeq::trecase

data(KLK1)
attach(KLK1)
names(KLK1)
str(KLK1[["Z"]]) # Just one gene and one SNP it seems in example 


# All matrices, rows = samples, cols = genes or SNPs
# Y = gene expression matrix
# Y1 = ASE gene expression matrix
# Y2 = ASE2 gene expression matrix
# X = covariate matrix: incl log reads per sample + PCs from log transformed expression matrix
# Z = phased genotype data: 0 AA, 1 AB, 3 BA, 4 BB

# eChr = chromosomes of gene expression traits - vector of integers 
# ePos = positions of Gex traits - vector of integers
# mChr = chromosomes of markers - vector of integers
# mPos = positions of markers - vector of integers

# Example has 65 participants, just one gene in Y/1/2, one genotype vector
# For eChr/pos just one numeric value entered as only one gene, thus chrom
# Not sure of the marker distance, do they mean end pos

# So in effect we can run per gene, with a genotype matrix per SNP as we have done before for CSeQTL

# TODO:

dim(SNP_test) # transform to samples by SNPs
dim(aseq_dat) # just subset, Y = aseq_dat$trec, Y1 = aseq_dat$hap1
dim(covars) # get relevant covars, re run on covar_modelling on residuals but no cell types perhaps 


# Set up -----------------------------------------------------------------------

library("stringr")
library("data.table")
library("dplyr")
library("smarter")
library("asSeq")

# Inherit arguments
# args <- commandArgs(trailingOnly = TRUE)
# cat("Inherited arguments:", paste(args, collapse = " "), "\n")
# array_file <- args[1]
# index <- as.integer(args[2])
# array <- read.csv(file = array_file)
# gene_id <- array[index, "ensembl_gene_id"]
# n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))

# For test
array <- read.csv(file = "/fh/working/hsu_l/Mari/cseqtl/scripts/eqtl_array_all_chr_rmdup.csv")
gene_id <- array[5, "ensembl_gene_id"]
n_cores <- 4

# Update array -----------------------------------------------------------------

# set up mart
library(biomaRt)
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
chroms <- (1:22)

# get coordinates
coord <- getBM(
  attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position', 'hgnc_symbol'),
  filters = c('chromosome_name', 'ensembl_gene_id'),
  values = list(chromosome_name=chroms, ensembl_gene_id=array$ensembl_gene_id),
  mart = ensembl)

coord <- inner_join(array, coord, by = c("ensembl_gene_id", "chromosome_name")) 

write.csv(coord,
          file = "/fh/working/hsu_l/Mari/cseqtl/scripts/trecase_array_all_chr_rmdup.csv", row.names = F, quote = F)

# Data Directories -------------------------------------------------------------

base_dir <- "/fh/working/hsu_l/Mari/cseqtl"

snp_dir <- file.path(base_dir, "genotype/phased/per_gene")
aseq_dir <- file.path(base_dir, "rnaseq/ASE/qc_genes/per_gene")
meta_dir <- file.path(base_dir, "metadata") 
out_dir <- file.path(base_dir, "eqtl/trecase_only")
setwd(out_dir)

# Load Covariates --------------------------------------------------------------

# Updated from 1_rnaseq_qc.R # See covar_modelling.R # update rna pcs
ids <- read.csv(file = file.path(meta_dir, "cseqtl_ids_rna_order.csv"))
# model_covars_base
covars <- read.csv(file = file.path(meta_dir, "cseqtl_sample_xx_vars.csv"))
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

# convert factor variables for matrix
model_covars$sct <- as.factor(model_covars$sct)
model_covars$ethnicity <- as.factor(model_covars$ethnicity)
model_covars$rna_batch <- as.factor(model_covars$rna_batch)
design <- model.matrix(~ ., data = model_covars) # remove intercept
design <- design[,-1]
dim(design)

# Load SNP data ----------------------------------------------------------------

# Genotype data made in aggregate, gathers SNPs per gene
setDTthreads(0L)
SNP <- fread(paste0(snp_dir,"/", gene_id, ".tsv"))
colnames(SNP) <- c("chrom", "pos", "ref", "alt", "genotype", "sample_id")

recode_genotype <- function(genotype) {
  ifelse(genotype == "0|0", 0,
         ifelse(genotype == "0|1", 1,
                ifelse(genotype == "1|0", 3,
                       ifelse(genotype == "1|1", 4, 5))))
} # NOTE: Different codes used compared with CSeQTL package

# Recode bases to CSeQTL codes
SNP[, `:=`(
  genotype_code = recode_genotype(genotype),
  snp_id = paste(chrom, pos, ref, alt, sep = ":")
)]

# Select the required columns, filter rows to include correct samples
SNP <- SNP[sample_id %in% ids$geno_id, .(sample_id, snp_id, genotype_code)]
SNP_bkp <- SNP

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
SNP <- t(SNP) # transform for aseq format


# AS-Seq data ------------------------------------------------------------------

# Load ASeq data 

# Aseq data made in Aseq format > aggregate_genes.R
aseq_file <- list.files(aseq_dir, pattern = gene_id)
aseq_dat <- fread(file.path(aseq_dir, aseq_file))
colnames(aseq_dat) <- c("sample_id", "trec", "hap1", "hap2", "hapN")
aseq_dat$ASREC <- aseq_dat$hap1 + aseq_dat$hap2

### Check orders ---------------------------------------------------------------

# Order aseq data, covars, and RHO matrix by SNP data
aseq_dat <- aseq_dat[match(row.names(SNP), aseq_dat$sample_id), ]
design <- design[match(row.names(SNP), row.names(design)), ]


# Gene loc file ----------------------------------------------------------------

# inf <- ENSG00000067064 10 1039152 1049119

chri <- 10 # inf[,2]
n_snps <- ncol(SNP_test)
n_genes <- 1 # ncol(aseq_dat$trec)


eChr = rep(chri, n_genes) # n_genes
eStart = 1039152 # inf[,3]
eEnd = 1049119 # inf[,4]
eExt = eEnd-eStart
pos = round((eStart+eEnd)/2)

mPos <- as.numeric(sapply(str_split(colnames(SNP_test), ":"), function(x) x[2])) # SNP/Marker positions
mChr = rep(chri, n_snps)


# Run Model --------------------------------------------------------------------

SNP_test <- SNP[ ,1:10]
# 10: 1,039,152-1,049,119

res.fil = "test2"
res.lon = sprintf("%s_eqtl.txt", res.fil)

trecase(Y = aseq_dat$trec,
        Y1 = aseq_dat$hap1,
        Y2 = aseq_dat$hap2,
        X = design,
        Z = SNP_test,
        eChr = eChr,
        mChr = mChr,
        ePos = pos,
        mPos = mPos,
        output.tag = res.fil,
        p.cut = 1)


# Output -----------------------------------------------------------------------

# res.fil = sprintf("v1_%s", inf[,1])



eqtl = read.table(res.lon, header=T, as.is=T)
eqtl[,1] = posi[eqtl[,2]]; colnames(eqtl)[1] = "Pos"
eqtl[,2] = vcfi[eqtl[,2]]
write.table(eqtl, res.lon, row.names=F, col.names=T, quote=F, sep="\t")
message("done version 1")

### Packages -------------------------------------------------------------------

# Bioconductor
library(BiocManager)
library(biomaRt)
library(Organism.dplyr)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(DOSE)
library(gwascat)
library(GenomicFeatures)
library(GenomicRanges)
library(LDlinkR)

# Data viz
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(viridis)
library(ggVennDiagram)

# Data cleaning
library(readODS)
library(forcats)
library(stringi)
library(stringr)
library(janitor)
library(data.table)
library(tidyr)
library(dplyr)
library(tidylog)

# Parallel
library(parallel)
library(doParallel)
library(foreach)

# options(timeout=1200)
# BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38", ask = FALSE, update = FALSE)
# library(SNPlocs.Hsapiens.dbSNP155.GRCh38)

# My eQTL data -----------------------------------------------------------------

# Directories
# base_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results"
# in_dir <- file.path(base_dir, "genotype/merged_study/eQTL/combined")
# plotdir <- file.path(base_dir, "plots/merged/eqtl") # scratch no longer available

in_dir="/fh/working/hsu_l/Mari"
setwd(in_dir)
load(file.path(in_dir, "cseqtl_results/mt_results_dec_thin.RData"))
in_dir="/fh/working/hsu_l/Mari"
eGenes$cell_type <- as.factor(eGenes$cell_type)


# Get annotation data ----------------------------------------------------------

# Convert SNP data to GRanges
snp_ranges <- GRanges(
  seqnames = eGenes$chr,
  ranges = IRanges(start = eGenes$pos, end = eGenes$pos),
  strand = rep("*", nrow(eGenes))
)

# Get gencode annotation
gtf_fn <- paste0(in_dir, "/SNEEP/resources/gencode.v43.annotation.gtf.gz")
exdb <- GenomicFeatures::makeTxDbFromGFF(file = gtf_fn, format = "gtf")
exons_list_per_gene <- GenomicFeatures::exonsBy(exdb, by = "gene")


### Filter for non-coding ------------------------------------------------------
# Directly check overlaps between SNPs and exons
coding_check <- overlapsAny(snp_ranges, unlist(exons_list_per_gene))
eGenes$genomic_context <- coding_check

non_coding <- filter(eGenes, genomic_context == FALSE) # 429,735 rows remaining

eGenes_nc <- filter(eGenes, snp %in% non_coding$snp)

# Assuming your main dataframe is named 'df'
df_list <- split(eGenes_nc, eGenes_nc$cell_type)

# Now, list_of_dfs is a list where each element is a dataframe.
# You can access each dataframe using the cell type as the key:
b_cells <- df_list$B_cells
cd4_t <- df_list$CD4_T_cells
cd8_t <- df_list$CD8_T_cells
mono <- df_list$Monocytes
neut <- df_list$Neutrophils
t_cells <- bind_rows(cd4_t, cd8_t)


# Create BED data frame with required columns (zero based start)
df <- mono

bed_data <- data.frame(
  chr = df$chr,
  start = df$pos-1,
  end = df$pos,
  var1 = df$ref,
  var2 = df$alt,
  rsid = "-",
  maf = -1
)

# Potentially duplicate SNPs, filter here
bed_data <- unique(bed_data)

write.table(bed_data, file = "SNEEP/snp_bed_data/mono_snps_eqtl_nc.bed",
            sep = "\t", quote = F, row.names = F, col.names = F)


# Get more annotation data -----------------------------------------------------

# List all feature types in the TxDb
featureTypes(exdb)

# Get a preview of genes
genes <- genes(exdb)
head(genes)

# Get a preview of exons
exons <- exons(exdb)
head(exons)

# Get a preview of transcripts
transcripts <- transcripts(exdb)
head(transcripts)

# GWAS catalogue/LD SNPs -------------------------------------------------------

ld_snps <- read.csv(file = "cseqtl_results/gwas_catalog_LD.csv")

## Format LD results for SNEEP -------------------------------------------------

# Update col names for BED format
ld_snps <- clean_names(ld_snps)

ld_snps <- ld_snps %>%
  separate(col = coord, into = c("chr", "pos"), sep = ":") %>%
  separate(col = alleles, into = c("ref", "alt"), sep = "/")

ld_snps$ref <- gsub("[()]", "", ld_snps$ref)
ld_snps$alt <- gsub("[()]", "", ld_snps$alt)
ld_snps$snp <- paste(ld_snps$chr, ld_snps$pos, ld_snps$ref, ld_snps$alt, sep = ":")
ld_snps$pos <- as.numeric(ld_snps$pos)

# Filter for eQTL-gwas SNP overlaps --------------------------------------------

eqtls <- unique(eGenes$snp)
eqtl_gwas <- filter(all_ld, snp %in% eqtls) # 14,442 SNPs
length(unique(eqtl_gwas$snp)) # 6612
length(unique(eqtl_gwas$rs_number)) 

# I believe there are duplicate rsIDs per position, diff minor allele freq in diff populations?
# And SNPs correlate with more than one seed:snp - remove dups, can refer back to all_ld later
eqtl_gwas <- eqtl_gwas %>% select(-seed_snp, -rs_seed, -distance, -dprime, -r2, -correlated_alleles, -maf) %>% distinct()

# one extra duplicate
eqtl_gwas <- eqtl_gwas[!duplicated(eqtl_gwas$rs_number), ]


# Get genomic context info -----------------------------------------------------

# Convert SNP data to GRanges
snp_ranges <- GRanges(
  seqnames = eqtl_gwas$chr,
  ranges = IRanges(start = eqtl_gwas$pos, end = eqtl_gwas$pos),
  strand = rep("*", nrow(eqtl_gwas))
)

# Get gencode annotation
gtf_fn <- paste0(in_dir, "/SNEEP/resources/gencode.v43.annotation.gtf.gz")
exdb <- GenomicFeatures::makeTxDbFromGFF(file = gtf_fn, format = "gtf")
exons_list_per_gene <- GenomicFeatures::exonsBy(exdb, by = "gene")

# Directly check overlaps between SNPs and exons
coding_check <- overlapsAny(snp_ranges, unlist(exons_list_per_gene))
eqtl_gwas$genomic_context <- coding_check

non_coding <- filter(eqtl_gwas, genomic_context == FALSE) # 5742 and  870 coding gwas_eQTLs



### Split by original cell type ------------------------------------------------

eGenes_gwas_nc <- filter(eGenes, snp %in% non_coding$snp)
# Add rs ids to eGenes
ld_rs <- eqtl_gwas %>% select(snp, rs_number)
eGenes_gwas_nc <- left_join(eGenes_gwas_nc, ld_rs, by = "snp")

# Assuming your main dataframe is named 'df'
df_list <- split(eGenes_gwas_nc, eGenes_gwas_nc$cell_type)

# Now, list_of_dfs is a list where each element is a dataframe.
# You can access each dataframe using the cell type as the key:
b_cells <- df_list$B_cells
cd4_t <- df_list$CD4_T_cells
cd8_t <- df_list$CD8_T_cells
mono <- df_list$Monocytes
neut <- df_list$Neutrophils
t_cells <- bind_rows(cd4_t, cd8_t)



df <- neut
# Create BED data frame with required columns (zero based start)
bed_data <- data.frame(
  chr = df$chr,
  start = df$pos-1,
  end = df$pos,
  var1 = df$ref,
  var2 = df$alt,
  rsid = df$rs_number
)

# Potentially duplicate SNPs, filter here
bed_data <- unique(bed_data)

write.table(bed_data, file = "SNEEP/snp_bed_data/neut_snps_eqtl_gwas_nc.bed",
            sep = "\t", quote = F, row.names = F, col.names = F)


# save.image(file = "gwas_eqtl_feb24.RData")



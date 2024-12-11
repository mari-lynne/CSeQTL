# GWAS catalogue comparisions
# Bioconductor
library(BiocManager)
library(biomaRt)
library(Organism.dplyr)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(DOSE)
library(gwascat)
library(GenomicRanges)

# Data viz
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(viridis)

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

# options(timeout=1200)
# BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38", ask = FALSE, update = FALSE)
# library(SNPlocs.Hsapiens.dbSNP155.GRCh38)

# My eQTL data -----------------------------------------------------------------

eqtls <- read.csv(file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT/mt_nov_results.csv")
# eqtls$seqnames <- gsub("chr", "", eqtls$chr)


# GWAS catalogue genes ---------------------------------------------------------

# Load GWAS catalog data
gwas_cat <- makeCurrentGwascat()
gwas_cat <- clean_names(as.data.frame(mcols(gwas_cat)))
# or
gwas_cat <- fread("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/resources/gwas_catalog_v1.0-associations_e113_r2024-10-21.tsv")
gwas_cat <- clean_names(gwas_cat)
gwas_cat <- gwas_cat %>% mutate(chr = str_c("chr", chr_id)) %>% rename(pos = chr_pos) %>% select(-chr_id)
gwas_cat$pos <- as.integer(gwas_cat$pos)

# Filter for relevant traits ---------------------------------------------------
# Just diseases
cvd_traits <- c("Stroke", "stroke", "artery", "atrial", "Atrial", "hypertension", "Hypertension", "blood prerssure", "Cardiovascular", "Myocardial", "Heart", "heart", "Thrombosis", "myocardial", "cardio", "Cardio", "Cardiac", "Blood pressure")

cvd_gwas <- gwas_cat[str_detect(gwas_cat$disease_trait, paste(cvd_traits, collapse = "|")), ]
cvd_gwas <- cvd_gwas %>% select(-date_added_to_catalog, -first_author, -date, -journal, -link, -replication_sample_size, -merged, -platform_snps_passing_qc, -cnv)

eqtl_cvd <- inner_join(cvd_gwas, eqtls, by = c("chr", "pos"))

# Split by cell type
eqtl_cvd_split <- split(eqtl_cvd, eqtl_cvd$cell_type)
list2env(eqtl_cvd_split, envir = .GlobalEnv)


# Epigenetic Data --------------------------------------------------------------

epi <- load(file ="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/snp_nexus/output/epi_summary.RData")

gwas_epi <- inner_join(eqtl_cvd, epi_features, by = c("chr", "snp"))
gwas_epi <- filter(gwas_epi, !str_detect(disease_trait, "kidney"))
table(gwas_epi$cell_type.x)
afr <- filter(gwas_epi, str_detect(initial_sample_size, "African"))

# Non neut gene assocs (with an epi feature)
gwas_epi2 <- filter(gwas_epi, cell_type.x != "Neutrophils")


# Interesting genes # With matching epigenomic features
# gwas_epi2 - C2, PPFIA4
# No matching epi feature (yet... but previously had them so need more roadmap data)
# gwas_epi SH3YL1

# take rows 3, 4 = CD8 T-cells, 8 = CD4 GSDMB example
interest <- gwas_epi2 %>% group_by(mapped_gene, cell_type.x) %>% slice_min(p_value)
cvd_genes <- c("C2", "PPFIA4", "ABO", "CXCL5", "IL18RAP", "GSDMB", "ZPBP2", "ARRB2", "FCHO1")
interest <- filter(interest, str_detect(mapped_gene, paste(cvd_genes, collapse = "|")))
interest <- interest %>% group_by(mapped_gene, cell_type.x) %>% slice_min(p_value, with_ties = FALSE)


# Neut genes
gwas_epi2 <- filter(gwas_epi, cell_type.x == "Neutrophils")
interest2 <- filter(gwas_epi2, context != "intron_variant" & context != "intergentic_variant")

afr <- filter(gwas_epi2, str_detect(initial_sample_size, "African"))

# SH3YL1 ORMDL3
cvd_genes <- c("PECAM1", "CXCL5")
interest2 <- filter(interest2, str_detect(mapped_gene, paste(cvd_genes, collapse = "|")))
interest2 <- filter(interest2, !str_detect(disease_trait, "kidney"))
interest2 <- interest2 %>% group_by(mapped_gene) %>% slice_min(p_bf, with_ties = FALSE)

final_genes <- as.data.frame(bind_rows(interest2, interest))
final_genes <- final_genes %>% select(cell_type.x, eta, p_bf, disease_trait, initial_sample_size, snps, context, reported_gene_s, mapped_gene, p_value, feature)
write.csv(final_genes, file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/cvd_genes_epi.csv", row.names = F)


afr <- afr %>% select(cell_type.x, eta, p_bf, disease_trait, initial_sample_size, snps, context, reported_gene_s, mapped_gene, p_value, feature)
write.csv(afr, file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/afr_genes_epi.csv", row.names = F)

# Search for CVD genes in other functional annotations -------------------------

# cvd_traits <- c("Stroke", "Ischemic stroke", "Large artery stroke", "Small vessel stroke", 
#                    "Myocardial infarction", "Atrial fibrillation", "Coronary artery disease", 
#                    "Coronary heart disease", "Heart failure", "Cardiac hypertrophy", "Blood pressure", 
#                    "Systolic blood pressure", "Diastolic blood pressure", 
#                    "Lipoprotein-associated phospholipase A2 activity and mass", 
#                    "Lipoprotein phospholipase A2 activity in cardiovascular disease", 
#                    "Blood protein levels in cardiovascular risk", "Blood pressure traits (multi-trait analysis)", 
#                    "Cholesterol", "HDL", "LDL", "Triglyceride", "cardiovascular", "lipoprotein", "stroke")
# cvd_traits <- str_to_lower(cvd_traits)
# pattern <- paste(cvd_traits, collapse="|")
# cvd_gwas <- gwas_cat[str_detect(gwas_cat$mapped_trait, pattern), ] # 72349 SNPs
# 
# filt_traits <- c("alcohol", "sleep", "migraine", "Alzheimer", "depression", "deprivation", "body fat percentage", "response")
# pattern <- paste(filt_traits, collapse="|")
# cvd_gwas <- cvd_gwas[!str_detect(cvd_gwas$mapped_trait, pattern), ] # 70058 SNPs
# 
# filt_traits <- c("smoking", "HIV", "drinkers", "Thiazide", "kidney")
# pattern <- paste(filt_traits, collapse="|")
# cvd_gwas <- cvd_gwas[!str_detect(cvd_gwas$disease_trait, pattern), ] # 68923

# GWAS catalogue comparisons
# Bioconductor
library(BiocManager)
library(biomaRt)
library(Organism.dplyr)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(DOSE)
library(gwascat)
library(GenomicRanges)
library(LDlinkR)

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

# Parallel
library(parallel)
library(doParallel)
library(foreach)

# options(timeout=1200)
# BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38", ask = FALSE, update = FALSE)
# library(SNPlocs.Hsapiens.dbSNP155.GRCh38)

# My eQTL data -----------------------------------------------------------------

# Directories
base_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results"
in_dir <- file.path(base_dir, "genotype/merged_study/eQTL/combined")
plotdir <- file.path(base_dir, "plots/merged/eqtl")

load(file.path(in_dir, "mt_results_dec_thin.RData"))

# GWAS catalogue genes ---------------------------------------------------------

# Load GWAS catalog data

gwas_cat <- fread("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/resources/gwas_catalog_v1.0.2-associations_e113_r2024-12-19.tsv")

# Tidy data 
gwas_cat <- clean_names(gwas_cat)
gwas_cat <- filter(gwas_cat, chr_id %in% as.character(c(1:22))) # Only include chr 1-22
gwas_cat <- gwas_cat %>% mutate(chr = str_c("chr", chr_id))
gwas_cat <- gwas_cat %>% dplyr::rename(pos = chr_pos) # %>% select(-chr_id)
gwas_cat$pos <- as.integer(gwas_cat$pos)

# Filter for relevant traits ---------------------------------------------------

# Use EFO terms

efo <- clean_names(fread("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/resources/gwas_catalog_trait-mappings_r2024-12-19.tsv"))

# Filter efo list for CVD related terms
efo <- filter(efo, parent_term == "Cardiovascular disease" | str_detect(efo_term, "lipids ratio"))

# Filter non relevant CVD diseases inc:
# autoimmune, arrhythmia/nervous system related disorders, trauma related outcomes, pre-eclampsia valve defects and congenital diseases.

efo_traits <- c("migraine", "Churg-Strauss syndrome", "Brugada syndrome",
                "familial sick sinus syndrome", "retinal vascular disorder", "skin vascular disease",
                "congenital anomaly of the great arteries",
                "familial long QT syndrome", "heart septal defect", 
                "Brugada syndrome", "hypotension", "arrhythmia", "Arrhythmia",
                "Chagas cardiomyopathy", "angioedema", "congenital", "Mitral", "fibromuscular dysplasia",
                "bacterial endocarditis", "atrial heart septal defect",
                "epistaxis", "diabetic foot", "raynaud disease", 
                "ectopy",
                "mitral valve disease", "aortic valve disease", "pulmonary valve disease",
                "tricuspid valve disease", "heart valve disease", "hemorrhoid", "Takotsubo cardiomyopathy",
                "conotruncal heart malformations", "elevation", "pericarditis",
                "CADASIL", "conduction system disorder", "congenital anomaly of cardiovascular system",
                "hypotension", "pericardial effusion", "congenital left-sided heart lesions",
                "chemotherapy", "sinoatrial node disorder",
                "sick sinus syndrome", "macrovascular complications of diabetes",
                "post-operative", 
                "transposition of the great arteries",
                "aortic valve insufficiency",
                "blood vessel injury",
                "bundle branch block",
                "cardiotoxicity",
                "Spontaneous coronary artery dissection",
                "cervical artery dissection",
                "coronary vasospasm",
                "endocarditis",
                "heart conduction disease",
                "inflammation of heart layer",
                "Pseudotumor cerebri",
                "rheumatic heart disease",
                "arteritis",
                "tachycardia",
                "ventricular outflow obstruction",
                "torsades de pointes",
                "vascular dementia",
                "chronic venous insufficiency",
                "fibrillation",
                "hemorrhage",
                "aortic coarctation",
                "retinal vein occlusion",
                "Arteritis",
                "myocarditis",
                "retinal vein occlusion",
                "premature cardiac contractions",
                "coronary restenosis",
                "pregnancy") 

efo <- efo[!str_detect(efo$efo_term, paste(efo_traits, collapse = "|")), ] # 701 traits
traits <- as.data.frame(unique(efo$efo_term)) 

gwas_efo <- inner_join(gwas_cat, efo, by = "disease_trait") %>% distinct()
traits <- as.data.frame(unique(gwas_efo$efo_term))


# Get GWAS genes
snp_genes <- unique(unlist(strsplit(gwas_efo$snp_gene_ids, ", ")))
gwas_genes <- unique(c(gwas_efo$upstream_gene_id, gwas_efo$downstream_gene_id, snp_genes)) # 4685 GWAS gene ids

gwas_eqtl_genes <- intersect(gwas_genes, eGenes$gene_name) # 1859 genes to finemap

# SUMMARY

# 8703 SNPs without lipid ratios
# 16,754 SNPs when including lipid ratios
# 56 vs 61 parent traits overall to examine

# With lipid ratios, 1859 GWAS-eGenes to examine


# For each set of GWAS SNPs, we have obtained correlated SNPs that are in linkage disequilibrium (LD) with any of the original SNPs. We used the LDProxy Tool [33] and extracted the proxy SNPs via their API functionality. Proxy SNPs from the European cohort with an R2 >  = 0.75 and within a window of ± 500 000 bp centered around the original SNP were added to the GWAS SNPs. The combined set of proxy and lead SNPs was used as input set to SNEEP.


# LD Proxy set up --------------------------------------------------------------

gwas_snps <- unique(gwas_efo$snps) # 7631 SNPs

# Do in batches of 1000 so not to overload API server
# Initialize a list to store batches
batch_list <- list()

# Calculate the number of batches needed
total_snps <- length(gwas_snps)
num_batches <- ceiling(total_snps / 1350) 

# Loop to create batches
for (i in 1:num_batches) {
  start_index <- (i - 1) * 1000 + 1
  end_index <- min(i * 1000, total_snps)
  batch_list[[i]] <- gwas_snps[start_index:end_index]
}

# Assign names to each batch for easier reference
names(batch_list) <- paste0("batch_", 1:num_batches)

# Remove first batch (already run as test)
batch_list <- batch_list[-1]

# Run LD Proxy VIP -------------------------------------------------------------

setwd("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/eQTL/combined/LDproxy")

# VIP API generated from https://ldlink.nih.gov/?tab=apiaccess and email authors for VIP access
# dbddeb68ba85

# Run in parallel

# Register parallel backend to use a single core per task
numCores <- 5 # detectCores()  # Or specify the number of cores you want to use
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Run batches in parallel without collecting results
foreach(i = seq_along(batch_list), .packages = 'LDlinkR') %dopar% {
  start <- Sys.time()  # Start time for performance monitoring
  
  LDproxy_batch(
    snp = batch_list[[i]],
    pop = c("EUR", "AFR"),
    r2d = "r2",
    token = "dbddeb68ba85",
    append = FALSE,
    genome_build = "grch38",
    win_size = "500000",
    api_root = "https://ldlink.nih.gov/LDlinkRest"
  )
  
  end <- Sys.time()  # End time for performance monitoring
  time_taken <- end - start  # Calculate time taken for the job
  
  cat("Completed batch", i, "in", time_taken, "seconds\n")
}

# Stop the parallel cluster after all jobs are done
stopCluster(cl)













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

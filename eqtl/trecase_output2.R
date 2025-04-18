
# Aims:

# Consolidate ASeq results files
# Correct for MT using eigenMT results

library("ggplot2")
library("ggpubr")
library("viridis")
library("stringr")
library("data.table")
library("dplyr")
library("tidylog")
library("janitor")

base_dir <- "/fh/working/hsu_l/Mari/cseqtl"

in_dir <- file.path(base_dir, "eqtl/trecase_only")
eigen_dir <- file.path(base_dir, "eqtl/eigenMT")
out_dir <- file.path(base_dir, "eqtl/ase_results")

setwd(in_dir)

## Load ASE results ------------------------------------------------------------

files <- list.files(path = in_dir, pattern = "eqtl.txt", full.names = TRUE)

# Initialize empty list to store data
eqtl_list <- list()

# Loop through files, read and append gene name
for (i in seq_along(files)) {
  file <- files[i]
  gene <- str_extract(file, "ENSG[0-9]{11}")
  df <- fread(file)
  df$gene_name <- gene
  eqtl_list[[i]] <- df
}

# Combine all into one data.table
all_eqtls <- rbindlist(eqtl_list, use.names = TRUE, fill = TRUE)

## Load Eigen MT ---------------------------------------------------------------

all_eigen <- clean_names(fread(file.path(eigen_dir, "All_snps_eigen.txt"))) %>%
  select(gene_name, tests) %>%
  rename(mEff = tests)

## Run MT correction -----------------------------------------------------------
all_eqtls <- rename(all_eqtls, p_nom = final_Pvalue)

calc_adjP <- function(eqtls, eigen) {
  # Add eigenMT number of tests (mEff)
  eqtls_mt <- inner_join(eqtls, eigen, by = "gene_name")
  
  # 1) Gene-level Bonferroni correction
  eqtls_mt <- eqtls_mt %>%
    mutate(p_bf = pmin(p_nom * mEff, 1))
  
  # 2) Genome-wide FDR correction using BH procedure on gene-level minimum p-values
  min_ps <- eqtls_mt %>%
    group_by(gene_name) %>%
    summarize(min_p_bf = min(p_bf), .groups = "drop") %>%
    arrange(min_p_bf) %>%
    mutate(rank = row_number(),
           n = n(),
           bh_thresh = 0.05 * (rank / n)) %>%
    filter(min_p_bf <= bh_thresh) %>%
    slice_max(rank, with_ties = FALSE)
  
  # Assign significance
  if (nrow(min_ps) > 0) {
    p_adj_cut <- min_ps$min_p_bf
    eqtls_mt <- eqtls_mt %>%
      mutate(is_eGene = p_bf < p_adj_cut)
  } else {
    eqtls_mt <- eqtls_mt %>%
      mutate(is_eGene = FALSE)
  }
  
  return(eqtls_mt)
}

all_eqtls <- calc_adjP(all_eqtls, all_eigen)

### Quick summaries ------------------------------------------------------------

eGenes <- filter(all_eqtls, is_eGene == TRUE)
length(unique(eGenes$gene_name)) # 8571
length(unique(eGenes$snp_id)) # 1,732,237 vs 416,052 CSeQTL

save.image(file = file.path(out_dir, "ase_results.RData"))


### Compare with CSeQTL --------------------------------------------------------
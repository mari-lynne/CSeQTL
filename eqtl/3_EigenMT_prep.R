# Aims:

# Update MT to match matrix qTL then prep for EigenMT

library(data.table)
library(stringr)
library(dplyr)
library(tidylog)
library(janitor)

library(ggplot2)
library(viridis)
library(patchwork)
library(Rmpfr)
library(qvalue)

# Directories ------------------------------------------------------------------

base_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results"
meta_dir <- file.path(base_dir, "metadata/merged")
workdir <- file.path(base_dir, "genotype/merged_study/eQTL/combined")
plotdir <- file.path(base_dir, "plots/merged/eqtl")

setwd(workdir)

load(file = "eqtl_fdr_results.RData")
write.csv(eqtls, file = file.path(workdir, "eqtls_all_cell_types.csv"), row.names = F, quote = F)
# Then clear dir and reload
eqtls <- read.csv(file = file.path(workdir, "eqtls_all_cell_types.csv"))

#eQTL results file (redo) ------------------------------------------------------

# Get final eQTL list - filter just for neuts (should still be all snps tested)

# Filter so just one SNP per cell type:
eqtls_final <- eqtls[!duplicated(eqtls$SNP), ]

head(eqtls_final[order(-eqtls_final$PVAL_final), ])


# Re calc zero p values
# 2) Recalculate 0 p values 
eqtls_final$PVAL_final <- pchisq(eqtls_final$LRT_eqtl, df = 1, lower.tail = FALSE)

# Extra SNPs
# extracted using generate snp_lists.sh as we didnt save ns p-vals (non-adj.p <0.05)

eqtls_final <- clean_names(eqtls_final)

# Split per chromosome
eqtls_split <- split(eqtls_final, eqtls_final$chr)

# Format for eigenMT
eqtls_split <- lapply(eqtls_split, function(df) {
  df %>%
    select(snp, gene_name, pval_final) %>%
    rename(`p-value` = pval_final)
})

write.csv(eqtls, file = file.path(workdir, "eqtls_all_cell_types.csv"), row.names = F, quote = F)
eqtls_final <- eqtls_final %>% 
  select(snp, gene_name, pval_final) %>%
  rename(`p-value` = pval_final)

fwrite(eqtls_final, file = file.path("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT", "eigen_all_qtls.txt"), sep = "\t", quote = FALSE, row.names = FALSE)


# Save results files 
results_dir = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased/snp_lists"

lapply(names(eqtls_split), function(name) {
  write.csv(eqtls_split[[name]], file = paste0(results_dir, "/", name, ".csv"), row.names = F)
})


# For NS file we have to get SNPs gene

# Then in array (so per chr)
#SBATCH array=1-22

chr=${SLURM_TASK_ARRAY_ID}
results_dir="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased/snp_lists"

cd ${results_dir}
sig_snp_file="chr${chr}.csv"
ns_snp_file="chr${chr}_snp_list.txt"

# Filter code...



# 3) Format for eigenMT
# 1. the SNP ID, 2. the gene ID, and 3. the cis-eQTL p-value. This file must contain a header line with the p-value column indicated by p-value.

eigen_in <- select(results_fdr, snp, gene_name, pval_final, cell_type) %>% rename(`p-value` = pval_final)

# Split per cell type
eigen_split <- split(eigen_in, eigen_in$cell_type)
eigen_split <- lapply(eigen_split, function(x) select(x, -cell_type))

# Save results files 
results_dir = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study"

lapply(names(eigen_split), function(name) {
  write.csv(eigen_split[[name]], file = paste0(results_dir, "/", name, ".csv"), row.names = F)
})

# Also make subdirs with same cell type names # TODO

### Filter array files ---------------------------------------------------------
array <- read.csv("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/eQTL/eqtl_array_all_chr_rmdup.csv")

neut_array <- filter(array, ensembl_gene_id %in% neut_snps$gene_name)
b_array <- filter(array, ensembl_gene_id %in% b_snps$gene_name)
mono_array <- filter(array, ensembl_gene_id %in% mono_snps$gene_name)
cd4_array <- filter(array, ensembl_gene_id %in% cd4_snps$gene_name)
cd8_array <- filter(array, ensembl_gene_id %in% cd8_snps$gene_name)

write.csv(neut_array, file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/eQTL/neut_eigen_array.csv", row.names = F)
write.csv(b_array, file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/eQTL/b_eigen_array.csv", row.names = F)
write.csv(mono_array, file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/eQTL/mono_eigen_array.csv", row.names = F)
write.csv(cd4_array, file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/eQTL/cd4_eigen_array.csv", row.names = F)
write.csv(cd8_array, file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/eQTL/cd8_eigen_array.csv", row.names = F)


# > cd4_array <- filter(array, ensembl_gene_id %in% cd4_snps$gene_name)
# filter: removed 11,068 rows (83%), 2,314 rows remaining
# > mono_array <- filter(array, ensembl_gene_id %in% mono_snps$gene_name)
# filter: removed 11,821 rows (88%), 1,561 rows remaining
# > cd8_array <- filter(array, ensembl_gene_id %in% cd8_snps$gene_name)
# filter: removed 10,215 rows (76%), 3,167 rows remaining
# > b_array <- filter(array, ensembl_gene_id %in% b_snps$gene_name)
# filter: removed 11,364 rows (85%), 2,018 rows remaining

# Delete extra files in CD4 array (from mistaken cd8 run)
cell_type <- "CD4_T_cells"
input_dir <- file.path(wd, cell_type)
files <- list.files(input_dir, pattern = "\\.txt$", full.names = FALSE)
genes <- str_extract(files, "ENSG\\d+")
extra <- files[!genes %in% cd4_array$ensembl_gene_id]
# 7 extra files
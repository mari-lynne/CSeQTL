
library(data.table)
library(dplyr)
library(tidylog)
library(stringr)
library(janitor)
library(smarter)
library(ggplot2)
library(viridis)
library(patchwork)
library(edgeR)

source("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/functions.R")

# Preprocess SCT RNA data ------------------------------------------------------

# dge_extra <- read_omic(name = "dge_black.rds", dir = in_dir)
# covar <- dge_extra$samples

# Define the file paths
main_dir <- c("/fh/scratch/delete90/kooperberg_c/mjohnson/sct-rna/data")
counts <- read_omic(name = "whi_topmed_to6_rnaseq_gene_reads.gct.gz", dir = main_dir)

# Filter duplicate genes
counts <- counts %>% distinct(Description, .keep_all = TRUE)
row.names(counts) <- counts$Description

genes <- counts[1:2]
counts <- counts[,-c(1:2)]

# Load Covariate data ----------------------------------------------------------

# Filtered for those with genotyping data in CSeQTL study
meta_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/metadata"

covar <- read.csv(file.path(meta_dir, "freeze12_pca.csv")) # geno outliers are filtered
covar$subject_id <- as.character(covar$subject_id)

# Filter covar to include only rows where 'subject_id' is in count ids
ids <- colnames(counts)
covar <- filter(covar, rnaseq_ids %in% ids) # 2063 samples
table(covar$ethnic) # multiethinc

# Ensure all ids are in covar
common_ids <- ids[ids %in% covar$rnaseq_ids]

# Filter the counts matrix to include only columns with common_ids
counts <- counts[, common_ids, drop = FALSE]

# Reorder 'covar' so that 'subject_id' matches the column order of 'counts'
covar <- covar[match(common_ids, covar$rnaseq_ids), ]
# Check if the reordering was successful
all.equal(covar$rnaseq_ids, colnames(counts))

counts <-as.matrix(counts)
dim(counts)

# Filter for previously qc passing genes (as in RNA data_merge)
dge_lls <- readRDS(file ="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/rnaseq/lls/lls_dgelist.rds")
dim(dge_lls)

# Sort gene info out
genes <- dge_lls$genes
genes <- genes %>% distinct(Description, .keep_all = TRUE) # FIlter dups

# Filter for matching genes
counts <- counts[row.names(counts) %in% genes$Description, ]
# rm(test)
dim(counts)
genes <- genes[row.names(counts) %in% genes$Description, ]
all.equal(genes$Description, row.names(counts))


##### Make DGE_list ------------------------------------------------------------

dge_sct <-  DGEList(counts=counts, samples=covar, genes=genes) 

# Filter spike genes
dge_sct <- dge_sct[!str_starts(row.names(dge_sct$genes), "ERCC"), ]

saveRDS(dge_sct, file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/rnaseq/sct/sct_cseqtl_dgelist.rds")

# Calc RNA pcs in rna_data_merge

# From script rna_data_merge, but should follow on directly after above anyways
# Set up =======================================================================

library(data.table)
library(dplyr)
library(tidylog)
library(stringr)
library(janitor)
library(smarter)
library(ggplot2)
library(viridis)
library(patchwork)
library(edgeR)

source("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/functions.R")

plot_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/plots/merged"
work_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/rnaseq" # merged_study/
meta_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/metadata"
out_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/rnaseq"

setwd(work_dir)

## Load DGE data ---------------------------------------------------------------

dge_sct <- readRDS(file = file.path(work_dir, "sct/sct_cseqtl_dgelist.rds")) # Updated in rnaseq_preprocess.R
dim(dge_sct)

dge_lls <- readRDS(file = file.path(work_dir, "lls/lls_dgelist.rds"))  
dim(dge_lls)


# Update IDs to common subject ids
colnames(dge_sct$counts) <- dge_sct$samples$subject_id
colnames(dge_lls$counts) <- dge_lls$samples$subject_id

# Merge RNA count matricies ----------------------------------------------------

dim(dge_sct$counts)

# Filter for matching rows
# Filter for overlapping genes, update sct to hgnc
head(row.names(dge_sct$counts)) # update to gene symbols
colnames(dge_sct$genes); colnames(dge_lls$genes)

# Only keep intersecting genes, also keeps them in the same order
dge_sct <- dge_sct[intersect(dge_sct$genes$Description, dge_lls$genes$Description), ]
dge_lls <- dge_lls[intersect(dge_sct$genes$Description, dge_lls$genes$Description), ]

# Extract count matricies, convert to df, left join
sct_counts <- as.data.frame(dge_sct$counts) %>% mutate(gene_symbol = row.names(dge_sct$counts))
lls_counts <- as.data.frame(dge_lls$counts) %>% mutate(gene_symbol = row.names(dge_lls$counts))

intersect(colnames(sct_counts), colnames(lls_counts)) # Dup samples
counts <- left_join(sct_counts, lls_counts, by = "gene_symbol")

# Filter dup counts (but poss use later to check batch differences)
extra <- select(counts, ends_with("x")) # removed lls ones
counts <- select(counts, !ends_with("x"))
rownames(counts) <- counts$gene_symbol
counts <- select(counts, -gene_symbol)

# Update with recent metatdata -------------------------------------------------

covar <- read.csv(file.path(meta_dir, "freeze12_pca.csv")) # geno outliers are filtered
covar$subject_id <- as.character(covar$subject_id)

# Filter covar to include only rows where 'subject_id' is in count ids
ids <- colnames(counts)
covar <- filter(covar, subject_id %in% ids) # 2063 samples

# Ensure all ids are in covar
common_ids <- ids[ids %in% covar$subject_id]

# Filter the counts matrix to include only columns with common_ids
counts <- counts[, common_ids, drop = FALSE]
# Reorder 'covar' so that 'subject_id' matches the column order of 'counts'
covar <- covar[match(common_ids, covar$subject_id), ]
# Check if the reordering was successful
all.equal(covar$subject_id, colnames(counts))

counts <- as.matrix(counts)

# Sort gene info out
genes <- dge_sct$genes

# Make DGE List ================================================================

dge <- DGEList(counts = counts, samples = covar, group = covar$rna_batch, genes = genes)
dim(dge)

# double check rna batch, maybe update :)

# Run PCA of data merged without BC ============================================

cpm <- cpm(dge)

pca_result <- prcomp(t(cpm), center=TRUE, scale=TRUE, rank=25)

# Basic plots
plot(pca_result)

# Add in variables of interest to cross examine pca by
pca_df <- as.data.frame(pca_result$x) %>%
  mutate(ethnic = as.factor(dge$samples$ethnicity),
         batch = as.factor(dge$samples$rna_batch),
         sct = as.factor(dge$samples$sct),
         age = dge$samples$age,
         bmi = dge$samples$bmi)

levels(pca_df$ethnic)
# pca_df$ethnic <- recode(pca_df$ethnic,
#                         '3' = 'A.A',
#                         '4' = 'His',
#                         '5' = 'E.A')

## Plot PCA --------------------------------------------------------------------

pc_pal <- c("#A3008D", "#53B4CC", "#FF0066")
batch_pal <- c("#3A3335", "#EB5E55")

# Batch
pca_df %>%
  ggplot(aes(PC1, PC2, color = batch)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = batch_pal) +
  theme_bw() +
  labs(color = "batch") +
  ggtitle("RNA-Seq PCA")
ggsave(filename = file.path(plot_dir, "rna_pc_bacthes_pc1_2.png"))

# Ethnicity
pca_df %>%
  ggplot(aes(PC1, PC2, color = ethnic)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = pc_pal) +
  theme_bw() +
  labs(color = "batch") +
  ggtitle("RNA-Seq PCA")

ggsave(filename = file.path(plot_dir, "rna_pc_ethnic_pc1_2.png"))

save.image(file = file.path(out_dir, "sct_lls_merged.RData"))
saveRDS(dge, file = file.path(out_dir, "sct_lls_merged.rds"))
saveRDS(dge, file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/sct_lls_merged.rds")

# PVE --------------------------------------------------------------------------

# Step 1: Calculate PVE from the prcomp object
pve_vals <- (pca_result$sdev^2) / sum(pca_result$sdev^2)

# Step 2: Create a data frame for PVE values
pve_df <- data.frame(
  PC = factor(seq_along(pve_vals)),
  V1 = pve_vals * 100  # Multiply by 100 to represent as a percentage
)

num_pcs = 20
pve_df <- pve_df[1:num_pcs,]

# PVE Plot
pve_pal <- viridis(
  n = num_pcs, 
  option = "plasma",    
  begin = 0.1,          
  end = 0.9,            
  direction = -1,
  alpha = 1
)

# Plot proportion of variance
ggplot(pve_df, aes(PC, V1, colour = PC)) +
  geom_line(group = "PC", colour = "#454b51") +
  geom_point(size = 2)  +
  ylab("\nProportion of Variance (%)\n") +
  scale_colour_manual(values = pve_pal) +
  theme_bw() + theme(legend.position = "none") +
  ggtitle("LLS + SCT RNAseq PCA")

ggsave(filename = file.path(plot_dir, "rna_merged_pve.png"))

# Update Covars ----------------------------------------------------------------

# extra cibersort performed locally, saved to CSeQTL one drive and copied over
covar2 <- read.csv(file.path(meta_dir, "sct_his_lls_covars.csv"))

# Add PCs
pca_df <- as.data.frame(pca_result$x)
colnames(pca_df) <- str_c(colnames(pca_df), "_rna")
pca_df$subject_id <- row.names(pca_df)
covar2$subject_id <- as.character(covar2$subject_id)
covar2 <- inner_join(covar2, pca_df, by = c("subject_id")) 

write.csv(covar2, file = file.path(meta_dir, "sct_his_lls_covars.csv"), row.names = F)

# Update DGElist ---------------------------------------------------------------
# Check if the reordering was successful
all.equal(covar2$subject_id, colnames(dge$counts))
dge <- DGEList(counts = dge$counts, samples = covar2, group = covar$rna_batch, genes = dge$genes)
dim(dge)

saveRDS(dge, file = file.path(out_dir, "sct_lls_merged.rds"))
saveRDS(dge, file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/sct_lls_merged.rds")
save.image(file = file.path(out_dir, "sct_lls_merged.RData"))



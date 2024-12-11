# Covariate checking:

# Aims:
# Run linear model on voom data with the required covars
# Run PCA on these residuals and save for CSeQTL

# Set Up -----------------------------------------------------------------------

## Packages --------------------------------------------------------------------

library(data.table)
library(stringr)
library(viridis)
library(ggplot2)
library(edgeR)
library(limma)
library(patchwork)
library(PCAForQTL)
library(ggpubr)
library(dplyr)
library(tidyr)
library(tidylog)
source("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/functions.R")

## Dirs and Vars ---------------------------------------------------------------

# Input Directories
base_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results"
rna_dir <- file.path(base_dir, "genotype/merged_study/rnaseq")
out_dir <- file.path(base_dir, "genotype/merged_study/rnaseq")
plot_dir <- file.path(base_dir, "plots/merged")
meta_dir <- file.path(base_dir, "metadata/merged")

setwd(rna_dir)

## Format data -----------------------------------------------------------------

dge <- readRDS(file = file.path(rna_dir, "sct_lls_merged.rds"))
dge$samples <- dge$samples %>% select(-"X.FID")
dge$samples <- dge$samples %>% select(-"norm.factors")

## Design Matrix ---------------------------------------------------------------

# Impute missing BMI values ~ 20 samples, use median
med <- median(dge$samples$bmi, na.rm = TRUE)
dge$samples$bmi <- if_else(is.na(dge$samples$bmi), med, dge$samples$bmi)
# 1 sample missing age data
med <- median(dge$samples$draw_age, na.rm = TRUE)
dge$samples$draw_age <- if_else(is.na(dge$samples$draw_age), med, as.numeric(dge$samples$draw_age))

# Check Factors!
dge$samples$rna_batch <- as.factor(dge$samples$rna_batch)
dge$samples$plate <- as.factor(dge$samples$plate)

# CSeQTL covariates (see variable importance check script)
covar_use <- dge$samples %>%
  select(ethnicity, sct, draw_age,
         Neutrophils, Monocytes_Macrophages, CD4_T_cells, CD8_T_cells, B_cells,
         PC1, PC2, PC3, PC4, PC5, rna_batch)

covar_use <- covar_use %>%
  rename(
    PC1_geno = PC1,
    PC2_geno = PC2,
    PC3_geno = PC3,
    PC4_geno = PC4,
    PC5_geno = PC5
  )

# Coefficients not estimable: NK_cells 

# Set categorical factors
covar_use$ethnicity <- as.factor(covar_use$ethnicity)
covar_use$sct <- as.factor(covar_use$sct)

design <- model.matrix(~ ., data = covar_use)
dim(design)

# Model data with limma voom --------------------------------------------------

dge <- calcNormFactors(dge)
v <- voom(dge, design, plot=TRUE)

fit <- lmFit(v, design)
efit <- eBayes(fit)
plotSA(efit, main="Final model: Mean-variance trend")


# Residual PCA ----------------------------------------------------------------
# Extract fitted values from limma gene expression model
fitted_values <- fitted.values(fit)

# Calculate residuals
# 'v' is an EList object with log2-CPM values
residuals <- v$E - fitted_values
residuals <- t(residuals)
dim(residuals)

# Data is already scaled, but double check if I should or not
res_pca <- prcomp(residuals, rank. = 25, center = TRUE, scale = TRUE)
dim(res_pca$x)

### Plot Residual PCA ----------------------------------------------------------

# Add in variables of interest to cross examine pca by
res_df <- as.data.frame(res_pca$x)
# colnames(res_df) <- str_c(colnames(res_df), "_res")
res_df$subject_id <- row.names(res_df)
plot_vars <- bind_cols(covar_use, dge$samples$subject_id)
plot_vars <- rename(plot_vars, subject_id = "...15")

res_df <- left_join(res_df, plot_vars, by = "subject_id")
# res_df <- res_df[order(res_df$ethnic, decreasing = TRUE), ]

res_df$ethnicity <- as.factor(res_df$ethnicity)
res_df$sct <- as.factor(res_df$sct)

plot_pca_grid(res_df,
              num_components = 7,
              var = "ethnicity",
              title = "RNASeq PCA of Residuals",
              save = TRUE,
              dir = plot_dir,
              file_name = "resid_pca_ethnic",
              width = 9, height = 6)

plot_pca_grid(res_df,
              num_components = 7,
              var = "Neutrophils",
              title = "RNASeq Residual PCA",
              save = TRUE,
              dir = plot_dir,
              file_name = "resid_pca_neut",
              width = 9, height = 6)

res_df <- res_df[order(res_df$sct), ]
plot_pca_grid(res_df,
              num_components = 7,
              var = "sct",
              title = "RNASeq Residual PCA",
              save = TRUE,
              dir = plot_dir,
              file_name = "resid_pca_sct",
              width = 9, height = 6)

# Bacth
plot_pca_grid(res_df,
              num_components = 7,
              var = "rna_batch",
              title = "RNASeq Residual PCA",
              save = TRUE,
              dir = plot_dir,
              file_name = "resid_pca_bacth",
              width = 9, height = 6)


### Plot PVE ------------------------------------------------------------------

# Step 1: Calculate PVE from the prcomp object
pve_values <- (res_pca$sdev^2) / sum(res_pca$sdev^2)

# Step 2: Create a data frame for PVE values
pve_df <- data.frame(
  PC = factor(seq_along(pve_values)),
  V1 = pve_values * 100  # Multiply by 100 to represent as a percentage
)

num_pcs = 25
pve_df <- pve_df[1:num_pcs,]
pve_df$cum_sum <- cumsum(pve_df$V1)

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
pve_plot <- ggplot(pve_df, aes(PC, V1, colour = PC)) +
  geom_line(group = "PC", colour = "#454b51") +
  geom_point(size = 2)  +
  ylab("\nProportion of Variance (%)\n") +
  scale_colour_manual(values = pve_pal) +
  theme_bw() + theme(legend.position = "none") +
  ggtitle("LLS Residual PCA")

pve_plot
str(pve_df)

pve_plot1|pve_plot
ggsave(filename = file.path(plot_dir, "res_rna_pve.png"))

# Investigate genes driving residual 1 variance 
loadings_pc1 <- res_pca$rotation[,1]
sort_loadings_pc1 <- sort(loadings_pc1, decreasing = TRUE)
head(sort_loadings_pc1)
tail(sort_loadings_pc1)
barplot(sort_loadings_pc1)
sort_loadings_pc1
rm(loadings_pc1)

# Save new covariate table -----------------------------------------------------

covar_new <- res_df[,-c(16:25)] # extra PCs
covar_new <- covar_new  %>% select(-subject_id)

# save.image(file = file.path(out_dir, "res_pca_modelling.RData"))


# Other residual plots ---------------------------------------------------------


# Plot residual versus fitted gene expression values
t_residuals <- t(residuals)

# Set seed
set.seed(123) 
n_genes <- 750
# Randomly sample indices for genes
gene_indices <- sample(nrow(t_residuals), n_genes)

# Subset residuals and fitted values matrices to the sampled genes
sampled_residuals <- t_residuals[gene_indices, ]
sampled_fitted <- fitted_values[gene_indices, ]

# Plot Residuals vs. Fitted Values for the sampled genes
plot(sampled_fitted, sampled_residuals,
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = paste("Fitted Gene Expression Values vs Residuals  (", n_genes, " genes)", sep = ""),
     pch = 16, col = rgb(0.1, 0.4, 0.6, 0.6))
abline(h = 0, col = "red")


# Check versus batch and other covariates

# Set seed
set.seed(123) 
n_genes <- 1250
# Randomly sample indices for genes
gene_indices <- sample(nrow(t_residuals), n_genes)
# Subset residuals and fitted values matrices to the sampled genes
sampled_residuals <- t_residuals[gene_indices, ]

var_residuals <- apply(sampled_residuals, 2, var)
covariate <- dge$samples$rna_batch

# Create a plot of residual variance by covariate levels
plot(covariate, var_residuals, 
     xlab = "RNA Batch",
     ylab = "Residual Variance",
     main = "Residual Variance by RNA Batch for Sampled Genes",
     pch = 16, col = rgb(0.1, 0.4, 0.6, 0.6))

t_test_result <- t.test(var_residuals ~ covariate) # Used Welch
# Display the t-test results
print(t_test_result)

# Recalculate residuals with PCs -----------------------------------------------

design2 <- model.matrix(~ ., data = covar_new)
dim(design2)

# Model data with limma voom 

v2 <- voom(dge, design, plot=TRUE)

fit2 <- lmFit(v2, design2)
efit2 <- eBayes(fit2)
plotSA(efit2, main="Final model: Mean-variance trend")

# Extract fitted values from limma gene expression model
fitted_values2 <- fitted.values(fit2)

# Calculate residuals
residuals2 <- v2$E - fitted_values2
dim(residuals2)

# Set seed
set.seed(123) 
n_genes <- 1250
# Randomly sample indices for genes
gene_indices <- sample(nrow(residuals2), n_genes)
# Subset residuals and fitted values matrices to the sampled genes
sampled_residuals <- residuals2[gene_indices, ]

var_residuals <- apply(sampled_residuals, 2, var)
covariate <- dge$samples$rna_batch

# Create a plot of residual variance by covariate levels
plot(covariate, var_residuals, 
     xlab = "RNA Batch",
     ylab = "Residual Variance",
     main = "Residual Variance by RNA Batch for Sampled Genes (post PC inclusion)",
     pch = 16, col = rgb(0.1, 0.4, 0.6, 0.6))

t_test_result <- t.test(var_residuals ~ covariate) # Used Welch
# Display the t-test results
print(t_test_result)

# Save covariates --------------------------------------------------------------

# Also ordering by genotyping data (cseqtl_test.R)

### Log library size -----------------------------------------------------------
# Each column is a sample, rows are genes, therefore colsums are the total lib size for that sample
dim(dge$counts)
lib.size <- log10(colSums(dge$counts))

# PCs and known covariates
covar_new <- res_df[,-c(16:25)] # minus extra PCs
covar_new <- covar_new [,-c(20:24)] # minus cell props
covar_new <- covar_new %>% select(-subject_id)
cseqtl_sample_xx_vars <- cbind(lib.size, covar_new)


write.csv(cseqtl_sample_xx_vars,
          file = file.path(meta_dir, "cseqtl_sample_xx_vars.csv"),
          row.names = F)

# Cell type proportions
rho <-  res_df %>%
  select(Neutrophils, Monocytes_Macrophages, CD4_T_cells, CD8_T_cells, B_cells)

write.csv(rho,
          file = file.path(meta_dir, "cseqtl_rho.csv"),
          row.names = F)

# IDs
ids <- dge$samples %>% select(subject_id, geno_id)
write.csv(ids,
          file = file.path(meta_dir, "cseqtl_ids_rna_order.csv"),
          row.names = F)



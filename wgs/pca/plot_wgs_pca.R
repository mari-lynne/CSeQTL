# Set up =======================================================================

# libraries
library(stringr)
library(data.table)
library(janitor)
library(dplyr)
library(tidylog)
library(viridis)

library(ggplot2) # 3.5.1
library(patchwork) # 1.1.1 > # devtools::install_github("thomasp85/patchwork")

# Directories
in_dir <- c("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased/pca")
plot_dir <- c("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/plots/merged/pca")
meta_dir <- c("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/metadata")

# Functions
'%!in%' <- function(x,y)!('%in%'(x,y))
source("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/functions.R")
source("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/pc_outlier_rm.R")

# Input data 
pca <- fread(file.path(in_dir, "lls_sct_concat_filtered_pca2_results.eigenvec")) # %>% dplyr::rename(IID = `#IID`)
eigenval <- fread(file.path(in_dir, "lls_sct_concat_filtered_pca2_results.eigenval")) # For PVE plot
covar <- read.csv(file.path(meta_dir, file = "LLS_SCT_metadata_jul24.csv")) 

# Filter meta data for just wgs samples ----------------------------------------

pca <- pca[pca$IID %in% covar$geno_id, ]

# check for missing and dup samples
miss <- pca[pca$IID %!in% covar$geno_id,]
dup <- pca[duplicated(pca$IID), ]

# Remove extra samples in covar
covar$IID <- covar$geno_id
covar <- covar[!duplicated(covar$IID), ] # removes extra lls entries 

pca_covar <- inner_join(pca, covar, by = c("IID"))

# Set factors for plot
pca_covar$ethnicity <- as.factor(recode(pca_covar$ethnic, '3' = 'A.A', '4' = 'His','5' = 'E.A'))
pca_covar$rna_batch <- as.factor(pca_covar$rna_batch)

# PVE --------------------------------------------------------------------------
pve <- data.frame(PC = 1:10, pve = (eigenval/ sum(eigenval) * 100))
pve$PC <- as.factor(pve$PC)
cumsum(pve$V1)

# Plot =========================================================================

wgs_plot_data <- pca_covar %>% select(starts_with("PC"), ethnicity, rna_batch) # ethnicity
wgs_plot_data <- as.data.frame(wgs_plot_data)

wgs_plot_data %>%
  ggplot(aes(x = PC8, y = PC9, color = ethnicity)) +
  geom_point(alpha = 0.8, size = 1.2) +
  theme_bw()

plot_pca_grid(wgs_plot_data,
              num_components = 7,
              var = "ethnicity",
              title = "PCA phased Freeze 12 WGS data (outliers removed)",
              save = TRUE,
              dir = plot_dir,
              file_name = "LLS_SCT_PCA_ethnic_outlier_rm",
              width = 9, height = 7)

# Remove outliers ==============================================================


# Identify outliers
# Get outliers without splitting
outliers <- get_outliers(wgs_plot_data)
# Get outliers with splitting by ethnicity
outliers_split <- get_outliers(wgs_plot_data, split = TRUE, var = "ethnicity", num_sd = 6)
outliers <- pca_covar[outliers_split, ]

# NWD163472, NWD190498, NWD869271, NWD683957, NWD635203, NWD277258, NWD874019, NWD398292, NWD863692, NWD750273

write.table(outliers$IID, file = file.path(meta_dir, "geno_outliers.txt"), col.names=F, row.names = F, quote = F)

# remove in plink and rerun pca

outliers$IID

# PVE --------------------------------------------------------------------------

pve$PC <- as.factor(pve$PC)

# Plot PVE (ancestry) ------------------------------------------------------------------
library(colorspace)
PC_Palette <-
  heat_hcl(
    10,
    h = c(30, -160),
    c = c(80, NA, 45),
    l = c(38, 79),
    power = c(0.85, 1.0),
    fixup = TRUE,
    gamma = NULL,
    alpha = 1
  )

ggplot(pve, aes(PC, V1, fill = PC)) + geom_point(stat = "identity") + ylab("Proportion of Variance (%)\n") + scale_color_manual(values = PC_Palette) + theme(legend.position = "none")

library(ggplot2)


ggplot(pve, aes(x = PC, y = V1, group = 1)) + 
  geom_line(color = "black") +
  geom_point(aes(fill = PC), color = "black", size = 3, shape = 21) +  # Line connecting the points
  ylab("Proportion of Variance (%)\n") + 
  scale_fill_manual(values = PC_Palette) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(filename = file.path(plot_dir, "pve.png"))



# Save data ====================================================================

write.csv(pca_covar, file = file.path(meta_dir, "freeze12_pca.csv"), row.names = F, quote = F)

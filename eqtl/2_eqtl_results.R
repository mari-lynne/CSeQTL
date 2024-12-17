# eQTL results

# Aims -------------------------------------------------------------------------

# Load and format all CSeQTL results
# Recalculate very low pvalues which were rounded to zero by R
# Basic Bonf. MT correction
# Some basic result summary tables and plots

# Packages ---------------------------------------------------------------------

# Run on "R version 4.2.0 (2022-04-22)"
library(data.table)
library(stringr)
library(dplyr)
library(tidylog)

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

# Clean Data -------------------------------------------------------------------

# Data is saved in text files per gene, with snps in rows, and cell type summary stats in columns (i.e wide format)
# Split data into results per cell type, as well as reconcatenating into one results DF (long format, but now grouped by cell type)
# In final results file, as well as concatenating cell type results, concatenate all gene results = one very long results file

### Check single file ----------------------------------------------------------
data <- fread(file = file.path(workdir, "ENSG00000171219_pvals.txt"))

# Define cell types
cell_types <- c("B_cells", "CD4_T_cells", "CD8_T_cells", "Neutrophils", "Monocytes_Macrophages")

# Initialize an empty list to store the subsets
split_dfs <- list()

# Loop through each cell type, subset, and filter the data in one step
for (i in cell_types) {
  # Get the relevant columns for the cell type
  matched_cols <- grep(i, names(data), value = TRUE)
  pval_col <- paste0("PVAL_final.", i)
  
  # Subset the data for the specific cell type
  split_dfs[[i]] <- data[, .SD, .SDcols = c("SNP", matched_cols)]
  
  # Calculate q-values using the Benjamini-Hochberg procedure
  split_dfs[[i]][, q_value := p.adjust(get(pval_col), method = "BH")]
}

# Access the B_cells data frame
bcells <- split_dfs[["B_cells"]]
cd8 <- split_dfs[["CD8_T_cells"]]
cd4 <- split_dfs[["CD4_T_cells"]]
neut <- split_dfs[["Neutrophils"]]
mono <- split_dfs[["Monocytes_Macrophages"]]

### Check all files ------------------------------------------------------------

files <- list.files(path = workdir, pattern = "pvals.txt", full.names = TRUE)
n_genes <- 13365

# Define the cell types of interest (has to match eqtl results df column suffixes)
target_cell_types <- c("B_cells", "CD4_T_cells", "CD8_T_cells", "Neutrophils", "Monocytes_Macrophages")

# Process each file using lapply
result_list <- lapply(files, function(i) {
  # Read the data from the file
  data <- fread(i)
  
  # Extract the gene name from the file name
  gene_name <- str_extract(i, "ENSG[0-9]{11}")
  message(paste("Processing eQTLs in gene", gene_name))  # Print the gene being processed
  
  # Add the gene name as a new column in the data table
  data[, Gene_Name := gene_name]
  
  # Initialize a list to store results for each cell type
  cell_type_results <- list()
  
  # Process each target cell type
  for (j in target_cell_types) {
    # Identify columns related to the current cell type
    matched_cols <- grep(j, names(data), value = TRUE)
    pval_col <- paste0("PVAL_final.", j)
    
    # Subset data for the specific cell type
    subset_data <- data[, .SD, .SDcols = c("SNP", "Gene_Name", matched_cols)]
    
    # Store the result for this cell type
    cell_type_results[[j]] <- subset_data
  }
  
  # Return the results for all cell types for this file
  return(cell_type_results)
})

# Initialize an empty list to store the combined results for each target cell type
cell_type_combined <- setNames(vector("list", length(target_cell_types)), target_cell_types)

# Combine results for each target cell type across all files
for (j in seq_along(target_cell_types)) {
  cell_type <- target_cell_types[j]
  
  # Extract and combine the filtered data for the current cell type from all files
  cell_type_combined[[cell_type]] <- rbindlist(lapply(result_list, function(res) res[[j]]))
}


bcells <- cell_type_combined[["B_cells"]]
cd8 <- cell_type_combined[["CD8_T_cells"]]
cd4 <- cell_type_combined[["CD4_T_cells"]]
neut <- cell_type_combined[["Neutrophils"]]
mono <- cell_type_combined[["Monocytes_Macrophages"]]

### Tidy separate cell type results --------------------------------------------

rm(filtered_data, cell_type_combined, data, df, result_list, split_dfs)

bcells$cell_type <- "B_cells"
cd8$cell_type <- "CD8_T_cells"
cd4$cell_type <- "CD4_T_cells"
neut$cell_type <- "Neutrophils"
mono$cell_type <- "Monocytes_Macrophages"

colnames(bcells) <- str_remove_all(colnames(bcells), ".B_cells")
colnames(cd8) <- str_remove_all(colnames(cd8), ".CD8_T_cells")
colnames(cd4) <- str_remove_all(colnames(cd4), ".CD4_T_cells")
colnames(neut) <- str_remove_all(colnames(neut), ".Neutrophils")
colnames(mono) <- str_remove_all(colnames(mono), ".Monocytes_Macrophages")

bcells$ETA <- as.numeric(bcells$ETA)
cd8$ETA <- as.numeric(cd8$ETA)
cd4$ETA <- as.numeric(cd4$ETA)
neut$ETA <- as.numeric(neut$ETA)
mono$ETA <- as.numeric(mono$ETA)

# Make final table  ------------------------------------------------------------

# Concat cell type results
eqtls <- bind_rows(bcells, cd8, cd4, neut, mono)
# Split SNP ID
eqtls <- eqtls[, c("CHR", "POS", "REF", "ALT") := tstrsplit(SNP, ":", fixed = TRUE)]
# Convert POS to numeric
eqtls <- eqtls[, POS := as.numeric(POS)]
length(unique(eqtls$SNP))
length(unique(eqtls$Gene_Name))

zero_ps <- filter(eqtls, PVAL_final == 0) # Handle separately with precision float numbers (46,000)
eqtls <- filter(eqtls, PVAL_final != 0)

# Some tidyting for MT script input
eqtls_all <- clean_names(bind_rows(zero_ps, eqtls))
eqtls_all <- eqtls_all %>% rename(p_nom = pval_final) %>% select(-adj_p)

library(readr)
write_csv(eqtls_all, file = "eqtls_all_snps.csv.gz")

## Adjust P values -------------------------------------------------------------

# Get number of SNPs per gene
# Recalculate p values using mpfr object to store high precision float numbers
# mpfr values don't work with functions - only simple calcs
# Therefore Bonferroni calculation has to be in base r

# Calculate adjusted p value per gene with simple Bonferroni (Note this is not the final MT method we used)
results_fdr <- eqtls %>%
  group_by(cell_type, Gene_Name) %>%
  mutate(PVAL_per_gene = pmin(PVAL_final * length(SNP), 1)) %>%   # max p value can be is 1
  filter(PVAL_per_gene < 0.05)


# Add column where SNPs are grouped by cell type then gene name and calculate the number of SNPs per gene
zero_ps <- zero_ps %>%
  group_by(cell_type, Gene_Name) %>%
  mutate(n_snps_per_gene = n()) # same as length(SNP)

summary(zero_ps$n_snps_per_gene)

zero_ps <- as.data.frame(zero_ps)

# Recalculate p-values that are rounded to zero
zero_ps$PVAL_prec <- pchisq(zero_ps$LRT_eqtl, df = 1, log.p = TRUE) # Get log p values
zero_ps$PVAL_prec <- mpfr(zero_ps$PVAL_prec, precBits = 128) # Save as mpfr for accuracy
zero_ps$PVAL_prec <- exp(zero_ps$PVAL_prec) # Convert p values to normal scale by exponent

# Multiply by n_snps for Bonferroni
zero_ps$PVAL_per_gene <- zero_ps$PVAL_prec * zero_ps$n_snps_per_gene

# Tidy dataframes to match
zero_ps$PVAL_per_gene <- formatMpfr(zero_ps$PVAL_per_gene, max.digits = 3)
PVAL_prec <- zero_ps$PVAL_prec # for back up calcs
zero_ps <- zero_ps %>% select(-adj.P, - n_snps_per_gene, -PVAL_prec)

results_fdr <- select(results_fdr, -adj.P)
results_fdr <- as.data.frame(results_fdr)
results_fdr$PVAL_per_gene <- format(results_fdr$PVAL_per_gene, scientific = TRUE, digits = 3)

results_fdr <- bind_rows(zero_ps, results_fdr)

# Summary Stats ----------------------------------------------------------------

# 444, 846 eqtls (also poss more due to extra cell type)
# 388, 872 unique SNPs

# Count the number of unique SNPs for each group
results_fdr_summary <- results_fdr %>%
  group_by(cell_type) %>%
  summarize(n_unique_eqtls = n_distinct(SNP))

results_fdr_summary <- results_fdr_summary %>%
  add_row(cell_type = "Total", n_unique_eqtls = sum(results_fdr_summary$n_unique_eqtls))

results_fdr_summary

# 265,143 rows remaining (1.5* as many CSeQTLs :))
table(results_fdr$cell_type)
table(results_fdr$eQTL) # trans-eqtls not calculated with trecase
table(results_fdr$CHR)

save.image(file = "eqtl_fdr_results.RData")

# Summary plots ----------------------------------------------------------------

# Filter for MT
bcells_fdr <- bcells %>% filter(adj.P < 0.05)
cd8_fdr <- cd8 %>% filter(adj.P < 0.05)
cd4_fdr <- cd4 %>% filter(adj.P < 0.05)
neut_fdr <- neut %>% filter(adj.P < 0.05)
mono_fdr <- mono %>% filter(adj.P < 0.05)

### Pie Chart ------------------------------------------------------------------
cell_counts <- c(
  `B-cells` = nrow(bcells_fdr),
  `CD8 T-cells` = nrow(cd8_fdr),
  `CD4 T-cells` = nrow(cd4_fdr),
  Neutrophils = nrow(neut_fdr),
  Monocytes = nrow(mono_fdr)
)

df <- data.frame(cell_type = names(cell_counts),
                 count = cell_counts)
# Calculate percentage
df$percentage <- df$count / sum(df$count) * 100

# Create the pie chart
ggplot(df, aes(x = "", y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "CS-eQTL Cell Type Distribution") +
  theme(legend.title = element_blank())

ggsave(filename = file.path(plot_dir, "pie_fdr.png"))

### Karyotpe plot --------------------------------------------------------------

library(karyoploteR)
library(rtracklayer)
library(GenomicRanges)

kp <- plotKaryotype(genome = "hg38")

# Convert SNP data to a GRanges object

# Convert SNP data to a GRanges object
snps_granges <- with(results, GRanges(seqnames = Rle(CHR),
                                   ranges = IRanges(start = POS, end = POS),
                                   SNP_ID = SNP,
                                   Gene_Name = Gene_Name,
                                   Cell_Type = cell_type,
                                   PVAL_adjust = adj.P))
table(snps_granges$Cell_Type)
chrs <- paste0("chr", 1:22)

# Basic Plot
kp <- plotKaryotype(genome="hg38", chromosome=chrs)
kpAddCytobandsAsLine(kp)
# kpPlotRegions(kp, data=snps_granges)

# Extract Cell_Type information from elementMetadata
cell_types <- mcols(snps_granges)$Cell_Type
# Filter snps_granges based on Cell_Type
b_cell_snps <- snps_granges[cell_types == "B_cells"]
t_cell_snps <- snps_granges[cell_types == "CD8_T_cells"|cell_types == "CD4_T_cells" ]
neut_snps <- snps_granges[cell_types == "Neutrophils"]
mono_snps <- snps_granges[cell_types == "Monocytes_Macrophages"]

# Plot SNPs for B cells using kpPlotRegions
kp <- plotKaryotype(genome="hg38", chromosome=chrs)
kpAddCytobandsAsLine(kp)
kpPlotRegions(kp, data = b_cell_snps, col = "#31688EFF", data.panel = 1, avoid.overlapping = TRUE)
kpPlotRegions(kp, data = t_cell_snps, col = "#35B779FF", data.panel = 1)
kpPlotRegions(kp, data = mono_snps, col = "#FDE725FF", data.panel = 1)
kpPlotRegions(kp, data = neut_snps, col = "#440154FF", data.panel = 1)
# Have to save in plot editor as pdf


# Demographic table ------------------------------------------------------------

covars <- read.csv(file = file.path(meta_dir, "cseqtl_sample_xx_vars.csv"))
covars <- read.csv(file = file.path(meta_dir, "sct_his_lls_covars.csv"))
covars$sct <- as.factor(covars$sct)
covars$rna_batch <- as.factor(covars$rna_batch)
covars$ethnicity <- as.factor(covars$ethnicity)

library(finalfit)

covars %>%
  summary_factorlist(dependent, explanatory, 
                     p=FALSE, add_dependent_label=FALSE,
                     column = TRUE) %>%
  add_col_totals() -> t1

knitr::kable(t1, align=c("l", "l", "r", "r", "r", "r"))


explanatory = c("ethnicity", "sct", "draw_age", "bmi")
dependent = 'rna_batch'

t1 <- covars %>%
  summary_factorlist(dependent, explanatory, 
                     p=FALSE, add_dependent_label=FALSE, 
                     column = TRUE, total_col = TRUE)

knitr::kable(t1, align=c("l", "l", "r", "r", "r", "r"))


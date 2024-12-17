

# Multiple Testing Correction - CSeQTL:

# 1a) EigenMT for Effective Number of Tests (meff): The EigenMT method calculates the effective number of tests (meff) for each gene by leveraging the correlation structure among single nucleotide polymorphisms (SNPs) within genes. Meff is calculated by identifying the smallest set of principal components (PCs) that capture the majority of the variance among the SNPs.

# 1b) Gene-level Bonferroni Correction: A Bonferroni correction is applied to each gene based on its meff. This yields a corrected p-value for each gene, accounting for the effective number of independent tests rather than the total number of SNPs.
                                        # P local = P value x Meff

# 2) Using the minimum p-value per gene calculate a gene level FDR threshold: Order all minimum p-values Determine sig threshold per ordered P-val.  - this also adjusts globally as considering n genes total

# 3) Identification of eGenes: Genes with at least one SNP locally adjusted p-value that meets the FDR significance threshold after Bonferroni correction are classified as "eGenes."

# 4) Extra SNPs are identified in Genes who's adjusted significance < gene level FDR

# Note: Do per cell type, as ultimately these are different tests with different distributions, also adjusting based on Neutrophils in B=cell SNPs where the significance is much higher in neut would create an unfairly stringent FDR threshold for B-cells

#### Packages ------------------------------------------------------------------

library(data.table)
library(dplyr)
library(stringr)
library(janitor)
library(tidylog)
library(readr)
library(ggplot2)
library(viridis)
library(ggVennDiagram)

# Data prep --------------------------------------------------------------------

# Directories
base_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results"
in_dir <- file.path(base_dir, "genotype/merged_study/eQTL/combined")
plotdir <- file.path(base_dir, "plots/merged/eqtl")

eqtls <- read_csv(file = file.path(in_dir, "eqtls_all_snps.csv.gz")) # Note this didn't save p-values of p_nom >= 0.05

# (Just look at qtls from cis/ASE model for now, add in others later)
eqtls <- filter(eqtls, e_qtl == "CIS") # removed 7% of tests!

# Recalculate 0 p values 
eqtls$p_nom <- pchisq(eqtls$lrt_eqtl, df = 1, lower.tail = FALSE)

# Check amount (now only 3640 zero values)
zero_ps <- filter(eqtls, p_nom == 0)
length(unique(zero_ps$gene_name))

## Load eigenMT results --------------------------------------------------------

wd <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT"
setwd(wd)

all_eigen <- clean_names(fread(file.path(wd, "All_snps_eigen.txt"))) %>% select(gene_name, tests)

# Calculate adjusted P values --------------------------------------------------

calc_adjP <- function(eqtls, eigen, Cell_Type = " ") {
  
  # Subset by cell type
  eqtls_cell <- eqtls %>% filter(cell_type == Cell_Type)
  # Add eigen MT number of tests (mEff)
  eqtls_cell_mt <- inner_join(eqtls_cell, eigen, by = "gene_name")
  
  # 1) Gene level Bonferroni correction ----------------------------------------
  eqtls_cell_mt <- eqtls_cell_mt %>%
    mutate(p_bf = p_nom * tests) %>%
    mutate(p_bf = pmin(p_bf, 1))  # Ensure p-values are capped at 1
  
  # 2) Genome-wide FDR correction ----------------------------------------------
  
  # Rank minimum snp p-values per gene
  min_ps <- eqtls_cell_mt %>%
    group_by(gene_name) %>%
    slice_min(p_bf, with_ties = FALSE) %>%
    ungroup()
  
  # Calculate FDR threshold across genes using BH procedure (alpha set to 0.05)
  min_ps <- min_ps %>%
    mutate(rank = rank(p_bf)) %>%
    mutate(bh_thresh = 0.05 * (rank / nrow(min_ps))) %>%
    select(gene_name, bh_thresh) 
  
  # Add  BH threshold for significant SNPs to DF
  eqtls_cell_mt <- left_join(eqtls_cell_mt, min_ps, by = "gene_name") 
  
  return(eqtls_cell_mt)
  
}

unique(eqtls$cell_type)
mono_mt <- calc_adjP(eqtls, all_eigen, Cell_Type = "Monocytes_Macrophages")

# 3) Identify eGenes -----------------------------------------------------------

all_mt <- bind_rows(b_mt, cd8_mt, cd4_mt, neut_mt, mono_mt) # combine cell type results

# Filter for significant SNPs by comparing Bonf. adjusted p values to the gene-level BH corrected threshold
eGenes <-  all_mt %>% filter(p_bf < bh_thresh)


save.image(file.path(in_dir, "mt_results_dec.RData"))


# Summarize data ---------------------------------------------------------------

gene_sum <- eGenes %>%
  group_by(cell_type) %>%
  summarise(n_genes = n_distinct(gene_name), .groups = 'drop')

# Pie chart isn't really a good way of summarizing as there are overlapping genes in each category
ggplot(gene_sum, aes(x = "", y = n_genes, fill = cell_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +  # Add white borders for clarity
  coord_polar(theta = "y") +  # Use polar coordinates
  theme_void() +  # Remove axis and grid
  labs(fill = "Cell Type", title = "\n          eGenes per Cell Type") + 
  theme(legend.position = "right",
        legend.text = element_text(size = 12)) +  # Increase legend text size
  scale_fill_viridis_d(direction = -1) +  # Apply reversed Viridis color scale
  geom_text(aes(label = n_genes), 
            position = position_stack(vjust = 0.5), 
            color = "white", 
            size = 5, 
            fontface = "bold")  # Add bold labels with n_genes values


### Venn diagram ---------------------------------------------------------------

# Plot n of overlapping eGenes between cell types

# Need to merge CD4 and CD8 results into T-cells for plotting purposes

eGenes2 <- eGenes %>% 
  mutate(cell_type = ifelse(cell_type %in% c("CD4_T_cells", "CD8_T_cells"), "T_cells", cell_type))


ven_genes <- list(
  `T cells` = eGenes2[eGenes2$cell_type == "T_cells", "gene_name"],
  `B cells` = eGenes2[eGenes2$cell_type == "B_cells", "gene_name"],
  Monocytes = eGenes2[eGenes2$cell_type == "Monocytes_Macrophages", "gene_name"],
  Neutrophils = eGenes2[eGenes2$cell_type == "Neutrophils", "gene_name"]
)

ggVennDiagram(ven_genes) + ggtitle("Cell-Type Specific eGenes\n") +
  ggplot2::scale_fill_gradient(low="violet",high = "gold")

ggsave(file = file.path(plotdir, "eGenes_venn.png"))




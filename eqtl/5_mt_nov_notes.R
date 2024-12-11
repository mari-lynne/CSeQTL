
# Check for missing genes 
eqtls_cell_mt <- left_join(eqtls, all_eigen, by = "gene_name")
no_eigen_gene <- filter(eqtls_cell_mt, is.na(tests))
length(unique(no_eigen_gene$gene_name))
# Summary per cell type
# B-Cells #22,000 SNPs with no matching gene in eigen results # 211 extra genes 
# Neutrophils # 84,417 SNPs, 213 missing genes
# TODO rerun EIGEN MT on 200+ missing genes
# For now, continue with MT on other genes, include eQTLs from this list if below genome wide significance maybe, and tack on at the end - BUT REDO 

# MT correction ================================================================

# MT adjustment per cell type then join at end
eqtls_cell <- filter(eqtls, cell_type == "Neutrophils")

# Step 1: Gene level Correction:
eqtls_cell_mt <- inner_join(eqtls_cell, all_eigen, by = "gene_name")

  # Calculate Bonferroni-corrected p-values for each SNP
  eqtls_cell_mt <- eqtls_cell_mt %>%
    mutate(p_bf = p_nom * tests) %>%
    mutate(p_bf = pmin(p_bf, 1))  # Ensure p-values are capped at 1

  # Get minimum p-value per gene after Bonferroni correction
  min_ps <- eqtls_cell_mt %>%
    group_by(gene_name) %>%
    slice_min(p_bf, with_ties = FALSE) %>%
    ungroup()

# Step 2: Genome wide FDR correction
  # Calculate FDR threshold across genes using BH procedure
  min_ps <- min_ps %>%
    mutate(rank = rank(p_bf)) %>%
    mutate(bh_thresh = 0.05 * (rank / nrow(min_ps))) %>%
    select(gene_name, bh_thresh)
 
# Step 3: Identify eGenes - genes with at least one SNP below threshold  
  # Filter for significant SNPs by comparing to gene-level BH threshold
  eqtls_cell_mt <- left_join(eqtls_cell_mt, min_ps, by = "gene_name") %>%
    filter(p_bf < bh_thresh)
  
b_qtls <- eqtls_cell_mt  
cd4_qtls <- eqtls_cell_mt
cd8_qtls <- eqtls_cell_mt
m_qtls <- eqtls_cell_mt


# Combine cell type results 
eqtls_cell_mt <- bind_rows(b_qtls, cd4_qtls, cd8_qtls, m_qtls, eqtls_cell_mt)
eGenes <- eqtls_cell_mt %>%
    group_by(gene_name) %>%
    slice_min(p_bf, with_ties = FALSE) 

# save.image(file = "mt_nov_results.RData")
# write.csv(eqtls_cell_mt, file = "mt_nov_results.csv", row.names = F, quote = F)

# Summary tables ===============================================================

# 1) Total n of SNPs sig after MT
# 2) Number of independent SNPs per gene (eigen results)
# 3) Number of eGenes per SNP

str(eqtls_cell_mt)

eqtls_cell_mt <- eqtls

eqtls_cell_mt %>%
  group_by(cell_type) %>%
  summarise(count = n(), .groups = 'drop')

gene_sum <- eGenes %>%
  group_by(cell_type) %>%
  summarise(n_genes = n_distinct(gene_name), .groups = 'drop')

gene_sum

library(ggplot2)

ggplot(gene_sum, aes(x = "", y = n_genes, fill = cell_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +  # Use polar coordinates
  theme_void() +  # Remove axis and grid
  labs(fill = "Cell Type", title = "Number of eGenes per Cell Type") + 
  theme(legend.position = "right") 


library(ggplot2)
library(viridis)

ggplot(gene_sum, aes(x = "", y = n_genes, fill = cell_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +  # Add white borders for clarity
  coord_polar(theta = "y") +  # Use polar coordinates
  theme_void() +  # Remove axis and grid
  labs(fill = "Cell Type", title = "Number of eGenes per Cell Type") + 
  theme(legend.position = "right") +
  scale_fill_viridis_d(direction = -1) +  # Apply Viridis color scale
  geom_text(aes(label = n_genes), 
            position = position_stack(vjust = 0.5), 
            color = "white", 
            size = 6,
            fontface = "bold")  # Add labels with n_genes values
library(ggplot2)
library(viridis)

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

ggsave(file = file.path(plot_dir, "eGenes_pie.jpg"))


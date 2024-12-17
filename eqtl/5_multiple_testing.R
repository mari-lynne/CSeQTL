

# Multiple Testing Correction - CSeQTL:

# 1) EigenMT for Effective Number of Tests (meff): The EigenMT method calculates the effective number of tests (meff) for each gene by leveraging the correlation structure among single nucleotide polymorphisms (SNPs) within genes. Meff is calculated by identifying the smallest set of principal components (PCs) that capture the majority of the variance among the SNPs.

# 2) Gene-level Bonferroni Correction: A Bonferroni correction is applied to each gene based on its meff. This yields a corrected p-value for each gene, accounting for the effective number of independent tests rather than the total number of SNPs.
                                        # P local = P value x Meff

# 3) Using the minimum p-value per gene calculate a gene level FDR threshold: Order all minimum p-values Determine sig threshold per ordered P-val.  - this also adjusts globally as considering n genes total

# 4) Identification of eGenes: Genes with at least one SNP locally adjusted p-value that meets the FDR significance threshold after Bonferroni correction are classified as "eGenes."

# 5) Extra SNPs are identified in Genes who's adjusted significance < gene level FDR

# Note: Do per cell type, as ultimately these are different tests with different distributions, also adjusting based on Neutrophils in B=cell SNPs where the significance is much higher in neut would create an unfairly stringent FDR threshold for B-cells

#### Packages ------------------------------------------------------------------

library(data.table)
library(dplyr)
library(stringr)
library(janitor)
library(tidylog)

# Data prep --------------------------------------------------------------------

# Using concatenated results from concat_eigen_results.sh - merge with OG results

# Directories
base_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results"
in_dir <- file.path(base_dir, "genotype/merged_study/eQTL/combined")
eqtls <- clean_names(read.csv(file = file.path(in_dir, "eqtls_all_snps.csv"), skipNul = T)) # Note this didn't save p-values of p_nom >= 0.05

eqtls <- eqtls

# Recalculate 0 p values 
eqtls$p_nom <- pchisq(eqtls$lrt_eqtl, df = 1, lower.tail = FALSE)

# Also just look at qtls from cis (ASE model for now, add in others later)
eqtls <- filter(eqtls, e_qtl == "CIS") # removed 7% of tests!

zero_ps <- filter(eqtls, p_nom == 0)
length(unique(zero_ps$gene_name))


### Load eigenMT results -------------------------------------------------------

wd <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT"
setwd(wd)

# Step 1) Gene-level Bonferroni Correction -------------------------------------

all_eigen <- clean_names(fread(file.path(wd, "All_snps_eigen.txt"))) %>% select(gene_name, tests)
# ns <- filter(all_eigen, !bf <0.05) # 18% 2,261 genes removed # 10,600 rows e-genes remaining

# Split into cell type lists

b_snps <- filter(eqtls, cell_type == "B_cells")

# Join each cell type df with all_eigen to get Meff

b_snps <- inner_join(b_snps, all_eigen, by = "gene_name") 

# For each cell type calculate the BH-adjusted P-values and take the minimum adjusted p-value per gene

b_snps <- b_snps %>%                                    
  mutate(p_bf = p_nom *tests) %>%
  mutate(p_bf = ifelse(p_bf > 1, 1, p_bf))
  
head(b_snps[order(b_snps$p_bf), ]) 

# Note: zero p-s i think this could skew fdr calcs, 
# Also should I at this stage filter for those without a BF < 0.05 min snp? - we are not testing those genes so why would we globally correct for them... leaving them in is more conservative...but actually we want the fdr stage to correct for all genes tested
# b_snps <- b_snps %>%
#   group_by(gene_name) %>%
#   filter(any(p_bf < 0.05))


# Step 2: Global FDR correction ------------------------------------------------


# Calculate FDR threshold per gene

# First get all the minimum p-values per gene
min_ps <- b_snps %>%
  group_by(gene_name) %>%
  slice_min(p_bf, with_ties = FALSE) 

# Next, order these p-values to get a per gene p-value threshold
min_ps <- min_ps %>% 
  mutate(rank = rank(p_bf)) %>%
  mutate(bh_thresh = 0.05 * (rank/nrow(min_ps))) %>% select(gene_name, bh_thresh)


# Step 3) Identify eGenes and extra SNPs passing threshold ---------------------

# Filter for BH p-values that are below bh_threshold per gene
b_fdr <- left_join(b_snps, min_ps, by = "gene_name") %>%
  group_by(gene_name) %>%
  filter(p_bf < bh_thresh)

# 2274 rows - 1069 genes 





# Get eGenes
all_eigen <- all_eigen %>% select(gene_name, tests) # %>% filter(bf <0.05) 

# Format results data 
all_snps <- inner_join(eqtls, all_eigen, by = "gene_name") %>% rename(p.value = pval_final) %>% select(-adj_p)

all_snps <- all_snps %>%
  mutate(p.adjust1 = p.value * tests,              # Bonferroni correction
         p.adjust1 = ifelse(p.adjust1 > 1, 1, p.adjust1))  # Cap values at 1

# Remove genes that dont have at least one SNP where p.adjust1 < 0.05
# Remove genes that don't have at least one SNP with p.adjust1 < 0.05
eGenes <- all_snps %>%
  group_by(cell_type, gene_name) %>%
  filter(any(p.adjust1 < 0.05)) %>% # Keep genes with at least one significant SNP
  ungroup()

length(unique(eGenes$gene_name)) # 12062 eGenes after round 1


# Step 2: Global FDR correction 

# Correct using the gene level adjusted p-values within each group (e.g., cell type)

all_snps <- all_snps %>% group_by(cell_type) %>%
  mutate(p.adjust2 = p.adjust(p.adjust1, method = "BH"))

sig_snps_fdr <- filter(all_snps, p.adjust2 < 0.05)
#  267,736 rows remaining sig SNPs remaining after FDR correction - # 8964 eGenes
length(unique(sig_snps_fdr$gene_name))   

# 3) Format data for SNP Nexus

# SNP Nexus --------------------------------------------------------------------

# chromosome	1	235781886	G	A	1

library(tidyr)

# Split the data by cell type into a list of data frames
snps_list <- split(sig_snps_fdr, sig_snps_fdr$cell_type)

# Define a function to process each subset, split into smaller chunks if needed, and save to text files
process_and_save <- function(data) {
  cell_type <- unique(data$cell_type)
  
  # Determine the number of rows in each chunk
  chunk_size <- 10000
  num_chunks <- ceiling(nrow(data) / chunk_size)
  
  # Split the data into smaller chunks if needed
  chunks <- split(data, ceiling(seq_len(nrow(data)) / chunk_size))
  
  # Process each chunk
  for (j in seq_along(chunks)) {
    formatted_data <- chunks[[j]] %>%
      separate(snp, into = c("chromosome_number", "position", "ref_allele", "alt_allele"), sep = ":", convert = TRUE) %>%
      mutate(chromosome = "chromosome", .before = 1,
             chromosome_number = str_extract(chromosome_number, "\\d+")) %>%
      select(chromosome, chromosome_number, position, ref_allele, alt_allele) %>%
      mutate(new_col = 1)
    
    # Save each chunk to a text file with a suffix indicating the chunk number
    file_name <- paste0(cell_type, "_", j, "_snp_nexus_in.txt")
    write.table(formatted_data, file_name, sep = "\t", quote = F, col.names = FALSE, row.names = FALSE)
  }
}

# Apply the function to each element in the list
lapply(snps_list, process_and_save)

# IN bash remove first column of all txt files
# for file in *.txt; do
# cut -f2- "$file" > temp && mv temp "$file"
# done

# Save data --------------------------------------------------------------------

write.csv(sig_snps_fdr, file = "sig_snps_fdr.csv", quote = F, row.names = F)

# Summarise new data -----------------------------------------------------------

plotdir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/plots/merged/eqtl"

table(sig_snps_fdr$cell_type)

# Filter to get top SNP per eGene, per cell type
filtered_results <- sig_snps_fdr %>%
  group_by(cell_type, gene_name) %>%
  filter(p.adjust2 == min(p.adjust2)) %>%
  ungroup() 

table(filtered_results$cell_type)

# Pie Chart
cell_type_df <- as.data.frame(table(filtered_results$cell_type))

# Rename the columns for clarity
names(cell_type_df) <- c("Cell_Type", "Count")

# Plotting the pie chart
pie(cell_type_df$Count, labels = cell_type_df$Cell_Type, 
    main = "Cell Type Distribution of eGenes", 
    col = rainbow(length(cell_type_df$Count))) 
ggsave(filename=file.path(plotdir,"cseqtl_eGenes_piechart.png"))


# Karyotype plot ---------------------------------------------------------------

library(karyoploteR)
library(rtracklayer)
library(GenomicRanges)

kp <- plotKaryotype(genome = "hg38")

# Convert SNP data to a GRanges object

# Convert SNP data to a GRanges object
snps_granges <- with(filtered_results, GRanges(seqnames = Rle(chr),
                                      ranges = IRanges(start = pos, end = pos),
                                      SNP_ID = snp,
                                      Gene_Name = gene_name,
                                      Cell_Type = cell_type,
                                      PVAL_adjust = p.adjust2))
table(snps_granges$Cell_Type)
chrs <- paste0("chr", 1:22)


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
kpPlotRegions(kp, data = b_cell_snps, col = "#31688EFF", data.panel = 1, avoid.overlapping = FALSE)
kpPlotRegions(kp, data = t_cell_snps, col = "#35B779FF", data.panel = 1, avoid.overlapping = TRUE)
kpPlotRegions(kp, data = mono_snps, col = "#FDE725FF", data.panel = 1, avoid.overlapping = TRUE)
kpPlotRegions(kp, data = neut_snps, col = "#440154FF", data.panel = 1, avoid.overlapping = TRUE)














# Checking meffs (differed between cell types)

# python /fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/resources/eigenMT/eigenMT.py --CHROM chr7 \
# --QTL "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT/temp/CD4_T_cells_ENSG00000002726_qtl_results.txt" \
# --GEN "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT//temp/ENSG00000002726_genotype_matrix.txt" \
# --GENPOS "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT/temp/ENSG00000002726_genotype_position.txt" \
# --PHEPOS "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT/temp/ENSG00000002726_phenotype_position.txt" \
# --cis_dist 500000 \
# --window 50 \
# --OUT "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT/ENSG00000002726_eigenMT_output_cd4_test1.txt"
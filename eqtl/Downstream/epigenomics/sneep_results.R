### Packages -------------------------------------------------------------------

# Bioconductor
library(BiocManager)
library(biomaRt)
library(Organism.dplyr)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(DOSE)
# library(gwascat)
library(GenomicFeatures)
library(GenomicRanges)
library(LDlinkR)

# Data viz
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(viridis)
library(ggVennDiagram)

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

# Load data ====================================================================

setwd("/fh/working/hsu_l/Mari")

# /fh/working/hsu_l/Mari/SNEEP/results/all_results

load(file = "cseqtl_results/gwas_eqtl_feb24.RData") 
load(file = "SNEEP/results/all_results/sneep_results.RData")


# Functional annotation of ALL CSeQTL loci -------------------------------------

eGenes <- filter(all_mt, is_eGene == TRUE)

# Convert SNP data to GRanges
snp_ranges <- GRanges(
  seqnames = eGenes$chr,
  ranges = IRanges(start = eGenes$pos, end = eGenes$pos),
  strand = rep("*", nrow(eGenes))
)

# Get gencode annotation
gtf_fn <- paste0(in_dir, "/SNEEP/resources/gencode.v43.annotation.gtf.gz")
exdb <- GenomicFeatures::makeTxDbFromGFF(file = gtf_fn, format = "gtf")
exons_list_per_gene <- GenomicFeatures::exonsBy(exdb, by = "gene")

# Directly check overlaps between SNPs and exons
coding_check <- overlapsAny(snp_ranges, unlist(exons_list_per_gene))
eGenes$genomic_context <- coding_check

non_coding <- filter(eGenes, genomic_context == FALSE) 

length(unique(non_coding$snp)) # 367,237 # 416,052 total unique
#88.3%


# SNEEP ------------------------------------------------------------------------

sneep <- fread("C:/Users/mjohnso5/Documents/cseqtl/sneep_results/neut_result.txt")

new_colnames <- c("SNP_position", "var1", "var2", "rsID", "MAF","peakPosition",
                  "TF", "strand", "TF-binding_position", "effectedPositionInMotif",
                  "pvalue_BindAff_var1", "pvalue_BindAff_var2", 
                  "log_pvalueBindAffVar1_pvalueBindAffVar2", "pvalue_DiffBindAff", "cell_type")

colnames(sneep) <- new_colnames


# MT adjustment
n_motifs <- length(unique(sneep$TF))
sneep$fdr_corrected_pvalue <- sneep$pvalue_DiffBindAff*n_motifs

b_sneep <- sneep
mono_sneep <- sneep
neut_sneep <- sneep
t_sneep <- sneep

b_sneep$cell_type <- "b_cells"
mono_sneep$cell_type <- "monocytes"
neut_sneep$cell_type <- "neutrophils"
t_sneep$cell_type <- "t_cells"

sig_sneep <- filter(all_sneep, fdr_corrected_pvalue < 0.05)
sig_sneep$SNP_TF_Pair <- paste(sig_sneep$SNP_position, sig_sneep$TF, sep = "_")

length(unique(sig_sneep$TF)) # 392 unique SNPs, 732 SNP-TF pairs, 257 TF's

top_sneep <- sig_sneep %>%
  group_by(cell_type, TF) %>%
  summarize(min_fdr = min(fdr_corrected_pvalue), .groups = 'drop') %>%
  group_by(cell_type) %>%
  arrange(min_fdr) %>%
  slice_head(n = 10) %>%
  ungroup()

top_sneep <- unique(top_sneep$TF)

sig_sneep <- sig_sneep %>%
  mutate(TF = ifelse(grepl("^TFAP2", TF), "TFAP2_Family", TF))

sig_sneep %>%
  filter(cell_type == "b_cells") %>%
  count(TF, sort = TRUE) %>%
  top_n(n = 10, wt = n) %>%
  ggplot(aes(x = "", y = n, fill = TF)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(fill = "Transcription Factor", title = "Top 10 rSNP Transcription Factors") +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))

sig_sneep %>%
  filter(cell_type == "monocytes") %>%
  group_by(TF) %>%
  summarize(min_fdr = min(fdr_corrected_pvalue), .groups = 'drop') %>%
  arrange(min_fdr) %>%
  slice_head(n = 10) %>%
  ggplot(aes(x = "", y = min_fdr, fill = TF)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(fill = "Transcription Factor", title = "Top 10 rSNP Transcription Factors in T-Cells") +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))

# save.image(file = "sneep_results.RData")

write.csv(sig_sneep, file = "cseqtl_results/summary_tables/sneep_results.csv", row.names = F)


# Load results -----------------------------------------------------------------

load(file = "ld_gwas_qtl.RData")

sig_sneep <- read.csv(file = "cseqtl_results/summary_tables/sneep_results.csv")
length(unique(sig_sneep$SNP_position)) # 392 SNEEP SNPs

sig_sneep <- select(sig_sneep, -rsID, -MAF)
sig_sneep <- sig_sneep %>%
  rename(TF.strand = TF.binding_position,
         TF.binding_position = strand,
         ref = var1,
         alt = var2) %>% select(SNP_position, ref, alt, cell_type, peakPosition, everything()) # TODO their ref/alt description is switched/wrong

# Link back to gene pair from eQTL_results

# Get snp pos (0 based so should be the second part of SNP_position)

sig_sneep <- sig_sneep %>% mutate(chr = str_extract(SNP_position, "chr[0-9XY]+"),
                                  pos = as.numeric(str_extract(SNP_position, "(?<=-)[0-9]+")),
                                  eQTL_snp = paste(chr, pos, ref, alt, sep = ":")) %>% select(-chr, -pos)

# Summarize functions 
eqtl_results <- eqtl_results %>%
  mutate(
    cell_type2 = case_when(
      cell_type == "B_cells" ~ "b_cells",
      cell_type == "CD4_T_cells" ~ "t_cells",  # Assuming you map both CD4 and CD8 T cells to "t_cells"
      cell_type == "CD8_T_cells" ~ "t_cells",
      cell_type == "Monocytes_Macrophages" ~ "monocytes",
      cell_type == "Neutrophils" ~ "neutrophils",
      TRUE ~ as.character(cell_type)  # Catch-all to keep other types unchanged if there are any
    )
  ) 

sig_sneep <- rename(sig_sneep, cell_type2 = cell_type)

sneep_genes <- left_join(sig_sneep, eqtl_results, by = c("eQTL_snp", "cell_type2"))

# Just GWAS loci ===============================================================

sneep_gwas <- filter(sneep_genes, cvd_gwas_qtl == TRUE)
unique(sneep_gwas$eQTL_snp)

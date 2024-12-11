# Set up -----------------------------------------------------------------------

# Data viz
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(viridis)

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

# Bioconductor
library(BiocManager)
library(biomaRt)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(DOSE)
library(gwascat)
library(GenomicRanges)

# base_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results"
# meta_dir <- file.path(base_dir, "metadata/merged")

base_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study"
workdir <- file.path(base_dir, "enrich")
plotdir <- file.path(base_dir, "plots/enrich")
outdir <- file.path(workdir)
setwd(workdir)

# Enrichment function
source("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/run_enrichment_function.R")

# Load data
# load(file = "eqtl_fdr_results.RData")
# rm(eqtls, neut, neut_fdr, mono, mono_fdr, bcells, bcells_fdr, cd4, cd4_fdr, cd8, cd8_fdr, PVAL_per_gene)

# Format meqtl results ---------------------------------------------------------------

meqtl <- meqtl %>%
  separate(col = mat_snp, into = c("CHR", "POS", "SNP_ID"), sep = ":")
meqtl$mat_gene <- sub(":.*", "", meqtl$mat_gene)

# check overlap in results
overlap <- results[results$SNP_ID %in% meqtl$SNP_ID, ]
length(unique(overlap$SNP_ID))

# check overlap in results
overlap <- results[results$Gene_Name %in% meqtl$gene, ]
length(unique(overlap$Gene_Name))
length(unique(meqtl$gene))
length(unique(results$Gene_Name))
meqtl <- filter(eqtls_full, !is.na(mat_gene))


# Enrichment function ==========================================================

library(dplyr)
library(tidylog)
library(Organism.dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(DOSE)

# Generate scr object beforehand
eqtls <- read.csv(file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT/mt_nov_results.csv")
genes <- unique(eqtls$gene_name) # 12199 genes

Sys.setenv(TMPDIR = "/home/mjohnso5/tmp")

src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")

# Options
# gene_col - name of column with gene info
# gene_type - type of gene id

df_list <- split(eqtls, eqtls$cell_type)

# Assign each dataframe in the list to a variable named after the cell type
list2env(df_list, envir = .GlobalEnv)

# write.table(Monocytes_Macrophages$Gene_Name, file = "mono_genes.txt")

enrich_output <- run_enrichment(data = Neutrophils,
                                gene_col = "gene_name",
                                gene_type = "ensembl",
                                name = "Neutrophil",
                                src_object = src)

# Unlist and assign objects to global environment:
for (ontology in names(enrich_output)) {
  assign(paste0(ontology), enrich_output[[ontology]])
}


# save.image(file = "enrich_results.RData")
# View(ego_BP_stroke@result)

# Plots ========================================================================

# Get all objects in the global environment that start with "ego"
enrich_objects <- ls(pattern = "^ego_MF")

# Loop through each ego object and plot
for (i in enrich_objects) {
  enrich_object <- get(i)
  if (any(enrich_object@result$p.adjust < 0.05)) {
    plt <- dotplot(enrich_object) + ggtitle(str_replace_all(i, "_", " "))
    print(plt)
    Sys.sleep(1.5)
    ggsave(filename = file.path(plotdir, paste0(i, ".png")))
  } else {
    cat("Skipping object", i, "as there are no significant results.\n")
  }
}

dotplot(ego_CC_meqtl)
ggsave(filename = file.path(plotdir, "ego_BP_meqtl.png"))
ggsave(filename = file.path(plotdir, "ego_MF_meqtl.png"))
ggsave(filename = file.path(plotdir, "ego_CC_meqtl.png"))

dotplot(ego_BP_Neutrophil, showCategory = 10) + ggtitle("Neutrophil")
ggsave(filename = file.path(plotdir, "bp_neut.png"))
dotplot(ego_BP_CD8_T_Cell, showCategory = 10) + ggtitle("T-cells")
ggsave(filename = file.path(plotdir, "bp_t_cell.png"))

cats <- c("viral process",
          "regulation of T cell activation",
          "T cell proliferation",
          "leukocyte cell-cell adhesion",
          "immune response-regulating signaling pathway",
          "leukocyte proliferation",
          "immune response-regulating cell surface receptor signaling pathway",
          "regulation of cell-cell adhesion",
          "viral life cycle",
          "regulation of leukocyte proliferation",
          "mononuclear cell proliferation")
dotplot(ego_BP_Monocyte, showCategory = cats) + ggtitle("Monocytes")
ggsave(filename = file.path(plotdir, "bp_monocyte.png"))


cats <- c("leukocyte cell-cell adhesion",
          "regulation of cell-cell adhesion",
          "myeloid cell differentiation",
          "response to virus",
          "immune response-regulating signaling pathway",
          "positive regulation of cytokine production",
          "positive regulation of leukocyte cell-cell adhesion",
          "regulation of T cell activation",
          "regulation of leukocyte cell-cell adhesion",
          "cytokine-mediated signaling pathway")

dotplot(ego_BP_B_Cell, showCategory = cats) + ggtitle("B Cells")
ggsave(filename = file.path(plotdir, "bp_bcell.png"))


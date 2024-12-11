# 2/10/23
# Aims
# Make gene coordinate file for splitting gene/snp matricies
# Updates - Jan 2024, include x chromosome - and only those with transcripts in both lls and sct
# Coordinates = 1 MB from 50000

# 1) Get gene coords ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(biomaRt)
library(dplyr)
library(tidylog)

# filter for genes with high/variable expression as examined in QC
base_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results"
meta_dir <- file.path(base_dir, "metadata")
study_genes <- read.csv(file = file.path(meta_dir, "gene_ids_outlierRM.csv"))
# Reformat ensembl id column (remove ./d)
study_genes$ensembl_gene_id <- str_replace(study_genes$Name, "\\..*", "")
chroms=c(1:22, "X")

# set up mart
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# get coordinates
coord <- getBM(
  attributes=c('chromosome_name','start_position','end_position','ensembl_gene_id', 'hgnc_symbol'),
  filters = c('chromosome_name', 'ensembl_gene_id'),
  values = list(chromosome_name=chroms, ensembl_gene_id=study_genes$ensembl_gene_id),
  mart = ensembl)

# 13877 genes - 13852 matched genes (there were no X chrom genes in study data)

# 2) Convert to eQTL coordinates (+-1MB or 500KB) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
coord$eqtl_start <- coord$start_position - 5e5
coord$eqtl_end <- coord$end_position + 5e5
# For now replace neg start positions (TODO check boundries)
coord$eqtl_start[coord$eqtl_start < 0] <- 1
# Order DF by chromosome (1:22, X) and then start position
coord$chromosome_name <- factor(coord$chromosome_name, levels = chroms)
coord <- coord[order(coord$chromosome_name, coord$start_position), ]

coord <- dplyr::select(coord, ensembl_gene_id, chromosome_name, eqtl_start, eqtl_end)

# Split into chromosome files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
chromosomes <- split(coord, coord$chromosome_name) # split into lists

for (chr in names(chromosomes)) {
  chr_df <- chromosomes[[chr]]
  file_name <- paste("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/resources/eqtl_500kb", chr, ".txt", sep="")
  write.table(chr_df, file=file_name, col.names=F, quote=F, sep="\t", row.names=F)
}

# without x
coord <- filter(coord, chromosome_name != "X")

# Save file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.table(coord, file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/resources/eqtl_coord.txt",
            col.names = F, quote = F, sep = "\t", row.names = F)

write.table(coord[5:25,], file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/resources/test_coord.txt",
            col.names = F, quote = F, sep = "\t", row.names = F)

# gene file names for agg_genes.sh
coord2 <- inner_join(coord, study_genes, by = "ensembl_gene_id")
write.table(coord2[, "Name"], file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/ASE/redo/qc_transcripts.txt",
            col.names = F, quote = F, sep = "\t", row.names = F)


# Make array for modelling script ----------------------------------------------

coord <- coord %>% select(ensembl_gene_id, chromosome_name)
chromosomes <- split(coord, coord$chromosome_name) # split into lists
chromosomes <- chromosomes[-23]


for (chr in names(chromosomes)) {
  chr_df <- chromosomes[[chr]]
  chr_df$index <- rep(1:nrow(chr_df))
  chr_df <- select(chr_df, index, ensembl_gene_id, chromosome_name)
  file_name <- paste("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/eQTL/eqtl_array_", chr, ".csv", sep="")
  write.csv(chr_df, file=file_name, quote=F, row.names=F)
}

# Slurm restart: ---------------------------------------------------------------

# see model_script_restart.sh for last gene file name extraction
# Remake array but minus the genes which have already been completed
# TODO need to check that some genes weren't cut off mid snp

last_genes <- fread("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/last_genes.txt", header = F)
coord2 <- coord[!coord$ensembl_gene_id %in% last_genes$V1, ]

chromosomes <- split(coord, coord$chromosome_name) # split into lists
chromosomes <- split(coord2, coord2$chromosome_name) # Restart (21 Feb)

for (chr in names(chromosomes)) {
  chr_df <- chromosomes[[chr]]
  chr_df$index <- rep(1:nrow(chr_df))
  chr_df <- dplyr::select(chr_df, index, ensembl_gene_id, chromosome_name)
  file_name <- paste("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/eQTL/eqtl_array_pt2_", chr, ".csv", sep="")
  write.csv(chr_df, file=file_name, quote=F, row.names=F)
}

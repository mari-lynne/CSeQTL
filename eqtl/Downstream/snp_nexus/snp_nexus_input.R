
library(data.table)
library(stringr)
library(janitor)
library(dplyr)
library(tidylog)

wd <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT"
setwd(wd)
load(file = "eqtl_multiple_testing.RData")


var <- "monocytes"
out_dir <- file.path("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/snp_nexus/output", var)


nexus <- as.data.frame(results$Monocytes_Macrophages)

nexus <-  nexus %>%
  select(chr, pos, ref, alt) %>%
  mutate(type = "Chromosome",
         strand = "1",
         chr = str_remove(chr, "chr")) %>%
  select(type, chr, pos, ref, alt, strand) # reorder


write.table(nexus, paste0(var,"_nexus.txt"), row.names = F, col.names = F, quote = F, sep = " ")


# FUMA Input ===================================================================

var <- "cd8_t_cells"
out_dir <- file.path("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/fuma")
setwd(out_dir)

# Liftover prep ----------------------------------------------------------------

# Use bed file (from epi_enrich.R), then download hg38 vcf file with RSIDs
# Then run bedtools intersect to match eQTL SNPs to the reference VCF - not very succesful
# liftover to hg37 instead

# Create a BED data frame with required columns
bed_data <- data.frame(
  chr = results_fdr$chr,
  start = results_fdr$pos,
  end = results_fdr$pos + 1,
  snp = results_fdr$snp
)
# Sort the data by chromosome and position
bed_data$chr <- factor(bed_data$chr, levels = paste0("chr", c(1:22)), ordered = TRUE)
bed_data <- bed_data[order(bed_data$chr, bed_data$start), ]

# eQTL SNP file:
write.table(bed_data, 
            file = file.path(out_dir, "eQTL_snps_38.bed"), 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)

#bedtools intersect -a eQTL_snps_38.bed -b /shared/biodata/ngs/GATK_hg38/dbsnp_146.hg38.vcf.gz -wa -wb > eqtl_rsids.txt


# Match with results -----------------------------------------------------------

snps_19 <- fread(file.path(getwd(), "lifted/eQTL_snps_19.bed"))
snps_19 <- distinct(snps_19)
names(snps_19) <- c("chr", "start", "end", "snp")
snps_19 <- select(snps_19, "start", "snp")

results_19 <- inner_join(results_fdr, snps_19, by = "snp")
fuma <- as.data.frame(filter(results_19, cell_type == "CD8_T_cells"))
# Select # fumma cols # use START for new POS
fuma <- select(fuma, chr, start, alt, ref, p_nom, eta)  # Use Beta as ETA, No SE calculated
colnames(fuma) <- c("CHROM","POS","A1","A2","P","ETA") # Note A1 = alt/affected allele

write.table(fuma, file.path(out_dir, "cd8_t_cells_fuma_input.txt"), row.names = F, quote = F, sep = "\t")

# get eigen results for lead snps

library(dplyr)

# Filter for the top SNP per gene_name, grouped by cell_type
top_snps <- results_19 %>%
  group_by(cell_type, gene_name) %>%
  arrange(p_nom) %>% # Sort by p-value within each group
  slice(1) %>%       # Select the top SNP (lowest p-value) in each group
  ungroup()          # Remove grouping

top_snps <- filter(top_snps, cell_type == "CD8_T_cells")
top_snps <- select(top_snps, chr, start, alt, ref, p_nom, eta)  # Use Beta as ETA, No SE calculated
colnames(top_snps) <- c("CHROM","POS","A1","A2","P","ETA") # Note A1 = alt/affected allele

write.table(top_snps, file.path(out_dir, "cd8_t_cells_lead_snps.txt"), row.names = F, quote = F, sep = "\t")

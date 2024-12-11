# SNP Nexus Results:

# Data viz
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(ggrepel)
library(ggridges)

# Data cleaning
library(readODS)
library(forcats)
library(stringi)
library(stringr)
library(janitor)
library(knitr)
library(data.table)
library(dplyr)
library(tidylog)
library(openxlsx)


cell_type <- "cd4_t_cells"

workdir <- file.path("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/snp_nexus/output", cell_type)
setwd(workdir)

plotdir <- c("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/plots")

# Load eQTL data
load(file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT/eqtl_multiple_testing.RData")


# Load nexus data --------------------------------------------------------------

# Download all results from SNPNexus, txt per annotation

file_names <- list.files(pattern = "*.txt") # Get file names
myfiles <- lapply(file_names, fread) #  read all files in temp (saved to list)
myfiles <- lapply(myfiles, clean_names) # clean column names
trim_names <- str_remove(file_names, ".txt") # tidy df names
names(myfiles) <- trim_names # Add names to list
list2env((myfiles), envir = .GlobalEnv) # Unlist files, save to global environment
KGen <- `1KGen`

# Do overall summaries, then A.A eQTL summaries

### Neutrophil data binding ----------------------------------------------------
# Read all files and save them in a list, also clean column names immediately
setDTthreads(0L)
file_names <- list.files(pattern = "\\.txt$") 
# remove mirna files, not much info
file_names <- file_names[!str_detect(file_names, "mir|regbuild|encode|pathway|ucs|ccd")]
myfiles <- lapply(file_names, fread)

# Extract prefixes (format is 'prefix_number.txt')
prefixes <- sapply(str_split(file_names, "_[0-9]+"), `[`, 1)

# Initialize an empty list to store combined data tables
combined_files <- list()

# Loop over unique prefixes to combine files with the same prefix
for (prefix in unique(prefixes)) {
  # Subset files that match the current prefix
  matched_files <- myfiles[prefixes == prefix]
  # Combine the matched files by binding rows
  combined_data <- rbindlist(matched_files)
  # Add the combined data table to the list, using the prefix as the name
  combined_files[[prefix]] <- combined_data
}

# Clean up names (remove ".txt" and numbers after underscores)
trim_names <- sapply(file_names, function(name) {
  parts <- unlist(str_split(name, "_[0-9]+\\.txt$"))
  parts[1]
})
names(combined_files) <- unique(trim_names)

# Export the list of data tables to the global environment
list2env(combined_files, envir = .GlobalEnv)
rm(myfiles, combined_files, matched_files)
KGen <- `1KGen`

# List all objects in the environment
all_objects <- ls()

# Apply clean_names to all data frames
for (obj_name in all_objects) {
  obj <- get(obj_name)
  if (is.data.frame(obj)) {
    assign(obj_name, clean_names(obj))
  }
}

save.image(file = paste0(cell_type, ".RData"))




# Epigenetics and genomic features ---------------------------------------------

# Filter for relevant cell type
# roadmap <- filter(roadmap, epigenome == "B cell")
# roadmap <- filter(roadmap, epigenome == "T cell")
# roadmap <- filter(roadmap, epigenome == "Neutrophil")
# roadmap <- filter(roadmap, epigenome == "Monocyte")

simple <- roadmap %>% select(variation_id, chromosome, feature_type_class)
simple <- distinct(simple) # 2096 rows

# Counting the number of occurrences for each feature_type_class
feature_counts <- simple %>%
  group_by(feature_type_class) %>%
  summarise(count = n()) %>%
  ungroup()

# Convert factor to use in pie chart (this will be used for the fill)
feature_counts$feature_type_class <- as.factor(feature_counts$feature_type_class)

# Create pie chart
ggplot(feature_counts, aes(x = "", y = count, fill = feature_type_class)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_viridis(discrete = TRUE, option = "plasma") + 
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = "white")) +
  labs(fill = "Feature Class", title = "Distribution of Epiginomic Feature Classes")

ggsave(filename = file.path(plotdir, paste0(cell_type, "epigenome_feature_classes.png")))

# Feature type histone markers
simple <- roadmap %>% select(variation_id, chromosome, feature_type)
simple <- filter(simple, str_starts(simple$feature_type, "H")) 

simple <- distinct(simple)
feature_counts <- simple %>%
  group_by(feature_type) %>%
  summarise(count = n()) %>%
  ungroup()

ggplot(feature_counts, aes(x = "", y = count, fill = feature_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_viridis(discrete = TRUE, option = "plasma") + 
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = "white")) +
  labs(fill = "Feature Type", title = "Distribution of Histone Markers")
ggsave(filename = file.path(plotdir, paste0(cell_type, "_histones_piechart.png")))

# Feature type - TF's
simple <- roadmap %>% select(variation_id, chromosome, feature_type)
simple <- filter(simple, !str_starts(simple$feature_type, "H"))

simple <- distinct(simple)
feature_counts <- simple %>%
  group_by(feature_type) %>%
  summarise(count = n()) %>%
  ungroup()

ggplot(feature_counts, aes(x = "", y = count, fill = feature_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_viridis(discrete = TRUE, option = "plasma") + 
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = "white")) +
  labs(fill = "Feature Type", title = "Distribution of TF's and other features")

ggsave(filename = file.path(plotdir, paste0(cell_type, "_tf_piechart.png")))
# Saving 7.26 x 4.52 in image


# 1KG data ---------------------------------------------------------------------
# Histogram of SNP frequencies x axis 0-1
# Replace "None" with 0, then make all the frequencies numeric variables
# I would then like to plot a histogram of SNP frequencies, per ancestry (0-1), with the ancestry histograms overlayed on top of eachother

# Convert freq columns to numeric
freq_cols <- c('afr_frequency', 'amr_frequency', 'eas_frequency', 'eur_frequency', 'sas_frequency')
KGen[, (freq_cols) := lapply(.SD, function(x) as.numeric(gsub("None", "0", x))), .SDcols = freq_cols]

# Box Plot
ggplot(KGen) +
  geom_density(aes(x = eas_frequency, color = 'East Asian'), alpha = 0.2, size = 0.9) +
  geom_density(aes(x = sas_frequency, color = 'South Asian'), alpha = 0.2, size = 0.9) +
  geom_density(aes(x = eur_frequency, color = 'European'), alpha = 0.2, size = 0.9) +
  geom_density(aes(x = afr_frequency, color = 'African'), alpha = 0.2, size = 0.9) +
  geom_density(aes(x = amr_frequency, color = 'American'), alpha = 0.2, size = 0.9) +
  labs(title = "eQTL Allele Frequencies by Ancestry", x = "Frequency", y = "Density", color = "Ancestry") +
  theme_bw()

ggsave(filename = file.path(plotdir, paste0(cell_type, "_1KG_ancestry_lines.png")))

# Ridgeline plot
long_KGen <- KGen %>%
  pivot_longer(cols = ends_with("frequency"),
               names_to = "Ancestry", values_to = "Frequency") %>%
  mutate(Ancestry = sub("_frequency", "", Ancestry))  # Tidy col names

# Assuming your data frame is named KGen
long_KGen <- KGen %>%
  pivot_longer(cols = ends_with("frequency"), names_to = "Ancestry", values_to = "Frequency") %>%
  mutate(Ancestry = sub("_frequency", "", Ancestry),  # Remove the '_frequency' part
         Ancestry = factor(Ancestry, levels = c("afr", "amr", "eas", "eur", "sas"),  # Ensure the order
                           labels = c("African", "American", "East Asian", "European", "South Asian")))  # Rename


ggplot(long_KGen, aes(x = Frequency, y = Ancestry, fill = Ancestry)) +
  geom_density_ridges() +
  scale_x_continuous(limits = c(0, 1)) +  # Assuming frequencies are between 0 and 1
  labs(title = paste(cell_type, " eQTL SNP Frequencies by Ancestry"), x = "Frequency", y = "Ancestry") +
  theme_minimal() +
  theme(legend.position = "none") 

ggsave(filename = file.path(plotdir, paste0(cell_type, "_1KG_ancestry_ridge.png")))

# GWAS associations ------------------------------------------------------------

# AFR ancestry GWAS
gwas_afr <- gwas[str_detect(gwas$population, "African"), ]
traits <- unique(gwas_afr$trait)
genes <- unique(gwas_afr$genes)

gwas_afr$p_value_numeric <- as.numeric(as.character(gwas_afr$p_value))  # Convert p-value to numeric
# Order by trait and p-value, then get the first row (which has the lowest p-value) for each trait
gwas_afr_ordered <- gwas_afr[order(trait, genes, p_value_numeric)]
# Collapse so only includes top SNP per trait
gwas_afr_collapsed <- gwas_afr_ordered[, .SD[1], by = trait]

gwas2 <- select(gwas_afr_collapsed, trait, variation_id, chromosome, allele_frequency, p_value, genes, population)
gwas2 <- select(gwas_afr_collapsed, trait, variation_id, chromosome, p_value, genes)
library(knitr)
gwas2$p_value <- formatC(gwas2$p_value, format = "e", digits = 2)
kable(gwas2, format = "markdown")

write.csv(gwas2, paste(cell_type, "gwas_afr_topsnp.csv"), row.names = F)

### Save/Load ------------------------------------------------------------------
# save.image(file = file.path(workdir, paste0(cell_type, "nexus_results.RData")))
# load(file = paste0(cell_type, "nexus_results.RData"))

# write list of genes for pathway enrichment
path_genes <- unique(ensembl$symbol) # 11,000 Genes

path_genes <- as.data.frame(path_genes)
write.csv(path_genes, file = paste0(cell_type, "gene_list.csv"), row.names = F)

# Merge results ----------------------------------------------------------------

clean_cols <- function(df) {
  # Tidy duplicate columns
  colnames(df) <- str_remove(colnames(df), ".x")
  # Remove duplicate column
  df <- df %>% dplyr::select(-matches("\\.y$"))
  return(df)
}

# Check for duplicate entries in ensemble and CADD data bases
# Some duplicates are because they have multiple entries per variant
# Therefore check to see if this is the case, if so concatenate entries where poss, otherwise just remove duplicate rows with stringr/distinct


# gen_coords - variant info, global freq, contigs
# KGEN ancestry specific frequencies
# Ensesmbl gene, variant data
# cadd, deletion scores

# cpg = cpg sites
# Roadmap = tissue specific epigenetic data

# GWAS - public GWAS results
# Clinvar - hereditary clinical associations

#### Simple duplicate removal --------------------------------------------------
library(stringi)
# Remove dups from other tables
gen_coords <- gen_coords[!stri_duplicated(gen_coords$variation_id), ]
cadd <- cadd[!stri_duplicated(cadd$variation_id), ]
cpg <- cpg[!stri_duplicated(cpg$variation_id), ]


#### Concatenate ---------------------------------------------------------------

# Coalesce duplicates from KGEN:
# Remove dup SNP with the lowest allele frequency. (Some dups just have 'none' recorded for freq)

KGen <- KGen %>%
  mutate(max_freq =
           pmax(afr_frequency, amr_frequency, eas_frequency, eur_frequency, sas_frequency, na.rm = TRUE)) %>%
  group_by(variation_id) %>%
  filter(max_freq == max(max_freq, na.rm = TRUE)) %>%
  ungroup()  %>%
  distinct()

# Coalesce predicted_function values for each variation_id
# First, remove transcripts info and duplicates
ensembl <- ensembl %>% dplyr::select(-transcript, -splice_distance) %>% distinct()
# Coalesce predicted function info per SNP
ensembl_2 <- ensembl %>%
  group_by(variation_id) %>%
  mutate(predicted_function = toString(unique(predicted_function))) %>% distinct()

### Merge results --------------------------------------------------------------

# Summary tables:

# 1) Gene and variant info:
# gen_coords, ensembl, 1KG, cadd, near_gens

# 2) Epigenetic:
# Roadmap + CPG

# 3) GWAS phenotype info:
# gwas + clinvar

# 1) Gene and variant info -----------------------------------------------------

# Select/Tidy Columns 
ensembl_2 <- select(ensembl_2,
                    variation_id, chromosome, position, variant, symbol, gene, predicted_function,
                    aa_position, aa_change, detail) %>% distinct()

KGen <- select(KGen, -db_snp)

# Merge tables
gene_var <- clean_cols(full_join(KGen, ensembl_2, by = "variation_id"))
gene_var <- clean_cols(full_join(gene_var, cadd, by = "variation_id"))
# entrez <- select(refseq, variation_id, entrez_gene)
# gene_var <- clean_cols(left_join(gene_var, entrez, by ="variation_id")) %>% distinct()
gene_var <- clean_cols(left_join(gene_var, sift, by = c("chromosome", "position")))  %>% distinct()

# Reorder
gene_var <- select(gene_var,
                       variation_id, chromosome, position, symbol, gene,
                       ref_allele, alt_allele, minor_allele, everything()) 


# 2) Epigenetics ---------------------------------------------------------------

# Filter for blood/immune cell expression
unique(roadmap$epigenome)

# Tidy tissue names
roadmap$epigenome <- str_remove_all(roadmap$epigenome, ("\\(PB\\)"))
roadmap$epigenome <- str_remove_all(roadmap$epigenome, ("Roadmap"))

# Filter for blood/immune cells
roadmap <-
  roadmap %>%
  filter(str_detect(roadmap$epigenome,
                    "B cell|CD4|CD8|CD34|T cell|T helper|Monocyte|neutrophil|
                    myeloid|Natural Killer|Spleen")==TRUE)

roadmap_bkup <- roadmap

# Tidy tissue names
roadmap <- roadmap %>%
  mutate(
    tissue_expression =
      case_when(
        grepl("common myeloid progenitor CD34 positive", epigenome) ~ "Common Myeloid Progenitor CD34+",
        grepl("CD4 positive alpha beta T cell", epigenome) ~ "CD4 T Cell",
        grepl("memory", epigenome) ~ "CD4 memory T-Cell",
        grepl("CD4 positive CD25 positive alpha beta regulatory T cell", epigenome) ~ "CD4 T-reg",
        grepl("naive thymus derived CD4 positive alpha beta T cell", epigenome) ~ "CD4 naieve T-Cell",
        grepl("B cell", epigenome) ~ "B Cell",
        grepl("T cell", epigenome) ~ "T Cell",
        grepl("Natural Killer cells", epigenome) ~ "NK Cell",
        grepl("neutrophil", epigenome) ~ "Neutrophil",
        grepl("T helper 17 cell", epigenome) ~ "Th-17 Cell",
        TRUE ~ epigenome
      )
  ) %>% select(-epigenome)

unique(roadmap$tissue_expression)
roadmap$tissue_expression <- str_trim(roadmap$tissue_expression)

# Update the start and end coordinates for each marker, grouped by variants and features type, across tissues
roadmap <- roadmap %>%
  group_by(variation_id, feature_type) %>%
mutate(region_start = min(region_start),
          region_end = max(region_end)) %>% distinct()

# Collapse duplicate histone markers by cell type (group by variation id))
roadmap <- roadmap %>%
  group_by(variation_id, feature_type) %>%
  mutate(tissue_expression = toString(unique(tissue_expression))) %>% 
  distinct()

# CPGs
epi <- clean_cols(full_join(roadmap, cpg, by = "variation_id"))
rm(roadmap_bkup)

# 3) GWAS ----------------------------------------------------------------------

gwas_2 <- clean_cols(full_join(gwas, clinvar, by = "variation_id")) %>% distinct()

# Maybe add a sheet for AFR GWAS too


# 4) eQTL results --------------------------------------------------------------

# # Matrix eQTL
# meqtl <- clean_names(fread("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/meqtl_results_corrected.txt"))
# colnames(meqtl) <- paste0("mat_", colnames(meqtl))
# meqtl$variation_id <- str_extract(meqtl$mat_snp, "[^:]*$")
# meqtl <- distinct(meqtl)
# 
# setDT(meqtl)  # Now meqtl is a data.table
# meqtl[, mat_fdr := as.numeric(mat_fdr)]
# meqtl <- meqtl[,.SD[which.max(mat_fdr)], by = variation_id]

# CSeQTL # Monocyte T_cells, B_Cells
# cqtl <- fread("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/Neutrophil_cseqtl_results.txt")

cqtl <- filter(results_fdr, cell_type == "CD4_T_cells") # Change
cqtl <- rename(cqtl, variation_id = snp,
               chromosome = chr,
               position = pos)
cqtl$chromosome <- paste0("chr", cqtl$chromosome)

# eqtl <- left_join(cqtl, meqtl, by = "variation_id") %>% distinct()

### Save xlsx ------------------------------------------------------------------

# Create a new workbook
wb <- createWorkbook()
# Add sheets to the workbook with your data frames
addWorksheet(wb, "eqtls")
addWorksheet(wb, "gene_variant_info")
addWorksheet(wb, "epigenetics")
addWorksheet(wb, "gwas_clinvar")
addWorksheet(wb, "gwas_afr")

writeData(wb, "eqtls", cqtl)
writeData(wb, "gene_variant_info", gene_var)
writeData(wb, "epigenetics", epi)
writeData(wb, "gwas_clinvar", gwas_2)
writeData(wb, "gwas_afr", gwas_afr)

# Save the workbook
saveWorkbook(wb, paste0(cell_type, "_eqtl_summary.xlsx"), overwrite = TRUE)

# rm(all_results, cadd, cell_sub, cell_sub2, check, clinvar, cpg, cqtl, ensembl, ensembl_2, entrez, KGen, meqtl_highest_fdr, refseq, results, test, gwas, roadmap)

# save.image(file = file.path(workdir, paste0(cell_type, "nexus_results_tidy.RData")))

# SQL --------------------------------------------------------------------------

library(RSQLite)
# Create SQLite database
db <- dbConnect(SQLite(), dbname = "snp_database.sqlite")
dbWriteTable(conn = db, name = "eqtl", value = cqtl, overwrite = TRUE)
dbWriteTable(conn = db, name = "gene_variant_info", value = gene_var, overwrite = TRUE)
dbWriteTable(conn = db, name = "epigenetics", value = epi, overwrite = TRUE)
dbWriteTable(conn = db, name = "gwas_clinvar", value = gwas_2, overwrite = TRUE)
dbWriteTable(conn = db, name = "gwas_afr", value = gwas_2, overwrite = TRUE)

# SQL Searching ----------------------------------------------------------------

# Search per table
snp_id <- 'rs111766855' 
query <- sprintf("SELECT * FROM eqtl WHERE variation_id = '%s'", snp_id)
results <- dbGetQuery(db, query)

# Search across tables
snp_id <- 'rs111766855' 
table_names <- dbListTables(db)
all_results <- list()

for (table in table_names) {
  query <- sprintf("SELECT '%s' as table_name, * FROM %s WHERE variation_id LIKE '%%%s%%'", table, table, snp_id)
  all_results[[table]] <- dbGetQuery(db, query)
}

# Initialize a variable for the joined data
joined_data <- NULL

# Iterate through the list of result tables
for (table in names(all_results)) {
  # If there's no joined data yet, set the first non-empty result as the starting point
  if (is.null(joined_data) && nrow(all_results[[table]]) > 0) {
    joined_data <- all_results[[table]]
  } else {
    # Join the data (this example uses a left join, change as needed)
    if (nrow(all_results[[table]]) > 0) {
      joined_data <- left_join(joined_data, all_results[[table]], by = "variation_id", suffix = c("", paste0("_", table)))
    }
  }
}

dbDisconnect(db)

# Combined data ================================================================

# workdir <- ("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/clean_results")
# setwd(workdir)
# load(file = "t_cellnexus_results.RData")
# load(file = "neut_nexus_results_tidy.RData")

b_qtl <- cqtl
b_epi <- epi
b_var <- gene_var
b_gwas <- gwas_2
b_afr <- gwas_afr

mono_cqtl <- cqtl
mono_epi <- epi
mono_var <- gene_var
mono_gwas <- gwas_2
mono_afr <- gwas_afr

cd4_qtl <- cqtl
cd4_epi <- epi
cd4_var <- gene_var
cd4_gwas <- gwas_2
cd4_afr <- gwas_afr

save.image("cseqtl_nexus_interim.RData")


# Merge Cell types =============================================================

gwas_afr_all <- bind_rows(b_afr, mono_afr, cd4_afr)




writeData(wb, "eqtls", eqtl)
writeData(wb, "gene_variant_info", gene_var)
writeData(wb, "epigenetics", epi)
writeData(wb, "gwas_clinvar", gwas_2)


# Create a new workbook
wb <- createWorkbook()
# Add sheets to the workbook with your data frames
addWorksheet(wb, "eqtls")
addWorksheet(wb, "gene_variant_info")
addWorksheet(wb, "epigenetics")
addWorksheet(wb, "gwas_clinvar")

writeData(wb, "eqtls", eqtl)
writeData(wb, "gene_variant_info", gene_var)
writeData(wb, "epigenetics", epi)
writeData(wb, "gwas_clinvar", gwas_2)

# Save the workbook
saveWorkbook(wb, paste0(cell_type, "_eqtl_summary.xlsx"), overwrite = TRUE)


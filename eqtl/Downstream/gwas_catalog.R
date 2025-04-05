# GWAS catalogue comparisons

### Packages -------------------------------------------------------------------

# Bioconductor
library(BiocManager)
library(biomaRt)
library(Organism.dplyr)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(DOSE)
library(gwascat)
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

clean_cols <- function(df) {
  # Tidy duplicate columns
  colnames(df) <- str_remove(colnames(df), ".x")
  # Remove duplicate column
  df <- df %>% dplyr::select(-matches("\\.y$"))
  return(df)
}

# My eQTL data -----------------------------------------------------------------

# Directories
# base_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results"
# in_dir <- file.path(base_dir, "genotype/merged_study/eQTL/combined")
# plotdir <- file.path(base_dir, "plots/merged/eqtl") # scratch no longer available

in_dir="/fh/working/hsu_l/Mari"
setwd(in_dir)
load(file.path(in_dir, "cseqtl_results/mt_results_dec_thin.RData"))
in_dir="/fh/working/hsu_l/Mari"
eGenes$cell_type <- as.factor(eGenes$cell_type)
# eGenes has been filtered for SNPs which pass MT correction at gene and genome wide level
# Has duplicates as SNPs counted in multiple cell types, and does include 'zero SNPs'

# GWAS catalogue genes ---------------------------------------------------------

# Load GWAS catalog data
gwas_cat <- fread(file.path(in_dir, "resources/gwas_catalog_v1.0.2-associations_e113_r2024-12-19.tsv"))

# Tidy data 
gwas_cat <- clean_names(gwas_cat)
gwas_cat <- filter(gwas_cat, chr_id %in% as.character(c(1:22))) # Only include chr 1-22
gwas_cat <- gwas_cat %>% mutate(chr = str_c("chr", chr_id))
gwas_cat <- gwas_cat %>% dplyr::rename(pos = chr_pos) # %>% select(-chr_id)
gwas_cat$pos <- as.integer(gwas_cat$pos)

# Filter for relevant traits ---------------------------------------------------

# Use EFO terms
efo_all <- clean_names(fread(file.path(in_dir, "resources/gwas_catalog_trait-mappings_r2024-12-19.tsv")))

# Filter efo list for CVD related terms
efo <- filter(efo_all, parent_term == "Cardiovascular disease"| disease_trait == "Asthma and cardiovascular disease") 
# | str_detect(efo_term, "lipids ratio"))

# Filter non relevant CVD diseases inc:
# autoimmune, arrhythmia/nervous system related disorders, trauma related outcomes, pre-eclampsia valve defects and congenital diseases.

efo_exclusions <- c("migraine", "Churg-Strauss syndrome", "Brugada syndrome",
                "familial sick sinus syndrome", "retinal vascular disorder", "skin vascular disease",
                "congenital anomaly of the great arteries",
                "familial long QT syndrome", "heart septal defect", 
                "Brugada syndrome", "hypotension", "arrhythmia", "Arrhythmia",
                "Chagas cardiomyopathy", "angioedema", "congenital", "Mitral", "fibromuscular dysplasia",
                "bacterial endocarditis", "atrial heart septal defect",
                "epistaxis", "diabetic foot", "raynaud disease", 
                "ectopy",
                "mitral valve disease", "aortic valve disease", "pulmonary valve disease",
                "tricuspid valve disease", "heart valve disease",
                "hemorrhoid", "Takotsubo cardiomyopathy",
                "conotruncal heart malformations", "elevation", "pericarditis",
                "CADASIL", "conduction system disorder", "congenital anomaly of cardiovascular system",
                "hypotension", "pericardial effusion", "congenital left-sided heart lesions",
                "chemotherapy", "sinoatrial node disorder",
                "sick sinus syndrome", "macrovascular complications of diabetes",
                "post-operative", 
                "transposition of the great arteries",
                "aortic valve insufficiency",
                "blood vessel injury",
                "bundle branch block",
                "cardiotoxicity",
                "Spontaneous coronary artery dissection",
                "cervical artery dissection",
                "coronary vasospasm",
                "endocarditis",
                "heart conduction disease",
                "inflammation of heart layer",
                "Pseudotumor cerebri",
                "rheumatic heart disease",
                "arteritis",
                "tachycardia",
                "ventricular outflow obstruction",
                "torsades de pointes",
                "vascular dementia",
                "chronic venous insufficiency",
                "fibrillation",
                "hemorrhage",
                "aortic coarctation",
                "retinal vein occlusion",
                "Arteritis",
                "myocarditis",
                "premature cardiac contractions",
                "coronary restenosis",
                "pregnancy",
                "peripartum cardiomyopathy") 

# Filter GWAS catalog data to only include entries which are in our CVD efo trait list

efo <- efo[!str_detect(efo$efo_term, paste(efo_exclusions, collapse = "|")), ] 
gwas_efo <- inner_join(gwas_cat, efo, by = "disease_trait") %>% distinct() 
exclude2 <- c("induced", "childhood", "Childhood") # Tidy more
gwas_efo <- gwas_efo[!str_detect(gwas_efo$disease_trait, paste(exclude2, collapse = "|")), ]
traits <- as.data.frame(unique(gwas_efo$disease_trait)) # 57 efo, 212 disease traits


### Filter to only look at GWAS/SNPs in eGenes ---------------------------------
# Get overlapping eQTL and GWAS genes
gwas_efo_wide <- gwas_efo %>%
  separate(snp_gene_ids, into = paste0("snp_gene_id", 1:3), sep = ",", remove = FALSE, fill = "right")

gwas_eqtl <-
  gwas_efo_wide %>% 
  filter(upstream_gene_id %in% eGenes$gene_name |
           downstream_gene_id %in% eGenes$gene_name |
           snp_gene_id1 %in% eGenes$gene_name |
           snp_gene_id2 %in% eGenes$gene_name |
           snp_gene_id3 %in% eGenes$gene_name)
length(unique(gwas_eqtl$snps))

gwas_genes <- unique(unlist(strsplit(gwas_efo$snp_gene_ids, ", ")))
gwas_genes <- unique(c(gwas_efo$upstream_gene_id, gwas_efo$downstream_gene_id, gwas_genes)) # 4088 GWAS gene ids
length(unique(gwas_efo$snps)); length(unique(gwas_genes))

gwas_eqtl_genes <- intersect(gwas_genes, eGenes$gene_name)
length(unique(gwas_eqtl$snps)); length(gwas_eqtl_genes)

gwas_efo <- gwas_eqtl # use this df from now on

# SUMMARY
# 2901 seed SNPs across 1610 CVD genes with eGenes in my data

# LD Proxy set up --------------------------------------------------------------

gwas_snps <- unique(gwas_efo$snps) # 8703 SNPs

gwas_snps <- miss # Rerun

# Do in batches of 1000 so not to overload API server
# Initialize a list to store batches
batch_list <- list()

# Calculate the number of batches needed
total_snps <- length(gwas_snps)
num_batches <- ceiling(total_snps / 1350) 

# Loop to create batches
for (i in 1:num_batches) {
  start_index <- (i - 1) * 1000 + 1
  end_index <- min(i * 1000, total_snps)
  batch_list[[i]] <- gwas_snps[start_index:end_index]
}

# Assign names to each batch for easier reference
names(batch_list) <- paste0("batch_", 1:num_batches)

# Remove first batch (already run as test)
batch_list <- batch_list[-1]

# Run LD Proxy VIP -------------------------------------------------------------

# setwd("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/eQTL/combined/LDproxy")

# VIP API generated from https://ldlink.nih.gov/?tab=apiaccess and email authors for VIP access
# dbddeb68ba85

# Run in parallel

# Register parallel backend to use a single core per task
numCores <- 2 # detectCores()  # Or specify the number of cores you want to use
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Run batches in parallel without collecting results
foreach(i = seq_along(batch_list), .packages = 'LDlinkR') %dopar% {
  start <- Sys.time()  # Start time for performance monitoring
  
  LDproxy_batch(
    snp = batch_list[[i]],
    pop = c("EUR", "AFR"),
    r2d = "r2",
    token = "dbddeb68ba85",
    append = FALSE,
    genome_build = "grch38",
    win_size = "500000",
    api_root = "https://ldlink.nih.gov/LDlinkRest"
  )
  
  end <- Sys.time()  # End time for performance monitoring
  time_taken <- end - start  # Calculate time taken for the job
  
  cat("Completed batch", i, "in", time_taken, "seconds\n")
}

# Stop the parallel cluster after all jobs are done
stopCluster(cl)

# LD Proxy Results -------------------------------------------------------------

# Run gwas_cat_file_cleaning.sh to filter for LD > 0.75

files <- list.files(path = file.path(in_dir, "ld_proxy/filtered"), pattern = "_filtered.txt$", full.names = TRUE) # 6461 SEED snp files (each with LD results) 

# Load SNP data and merge into data frame
data_list <- lapply(files, function(file) {
  fread(file)}) 

ld_snps <- rbindlist(data_list)

# Tidy LD_SNP data frame
ld_snps <- clean_names(ld_snps)
ld_snps <- ld_snps %>%
  rename(rsid_ld = rs_number, coord_ld = coord) %>%
  select(-maf, -distance, -dprime, -r2, -correlated_alleles) %>%
  distinct() # some snps corr with multiple seed snps, prevents dups later. # TODO add LD/distance info back in

write.csv(ld_snps, file = "cseqtl_results/gwas_catalog_LD.csv", row.names = F)

# Missing LD SNPs --------------------------------------------------------------

# Check for missing snps rerun LD proxy on them (note reran - these were ones with no match in 1000G db or monoallelic)
# Add missing snps to ld_snp list
seed_snps <- unique(ld_snps$seed_snp)
miss <- gwas_efo %>% filter(!snps %in% seed_snps) # SNPs not in LD run but in original GWAS list 
miss <- unique(miss$snps) # 838

# just rerun the asthma cvd ones
extra <- filter(gwas_eqtl, efo_term == "asthma")
miss <- intersect(miss, extra$snps) # 99 :)

setwd("/fh/working/hsu_l/Mari/ld_proxy/extra_snps")

start <- Sys.time()  # Start time for performance monitoring

LDproxy_batch(
  snp = miss,
  pop = c("EUR", "AFR"),
  r2d = "r2",
  token = "dbddeb68ba85",
  append = FALSE,
  genome_build = "grch38",
  win_size = "500000",
  api_root = "https://ldlink.nih.gov/LDlinkRest"
)

end <- Sys.time()  # End time for performance monitoring
time_taken <- end - start  # Calculate time taken for the job

cat("Completed batch", i, "in", time_taken, "seconds\n")

# Add to LD snps

# Get allele info from dbSNPs_sorted.txt 
extra <- gwas_efo %>% filter(!snps %in% seed_snps) %>% select(snps)
write.table(extra, file = file.path(in_dir, "SNEEP/resources/extra_snps.txt"),
            row.names = F, quote = F, col.names = F, sep = "\t")

# ml parallel
# parallel --pipepart -a dbSNPs_sorted.txt -k grep -Fwf extra_snps.txt > extra_snps_coord.txt

extra_coord <- fread("SNEEP/resources/extra_snps_coord.txt") # 559
# Format to match ld_snps
extra_coord <- extra_coord %>%
  rename(rs_number=V7, chr=V2, pos=V4, ref=V5, alt=V6) %>%
  select(rs_number,chr,pos,ref,alt) %>%
  mutate(maf=NA, distance = 1, dprime = 1, r2 = 1,
         correlated_alleles = NA, forg_edb = NA, regulome_db = NA, `function` = NA, seed_snp = rs_number,
         snp=paste(chr, pos, ref, alt, sep = ":"))
  
ld_snps <- bind_rows(ld_snps, extra_coord)
ld_snps$rs_seed <- paste(ld_snps$rs_number, ld_snps$seed_snp, sep = ":")
all_ld <- ld_snps


### Join LD + GWAS table -------------------------------------------------------

# Join back LD results to original GWAS loci and thus trait using seed snp
# Should give multiple entries per GWAS trait with corresponding LD snps
# Then join back to CSeQTL loci via LD snp coord

# Match up LD snps with the CVD gwas catalog data by SEED SNP rs ID
gwas_efo <- gwas_efo %>%
  select(-snp_id_current, -chr_id) %>%
  rename(seed_chr = chr, seed_pos = pos, seed_snp = snps) %>%
  mutate(seed_coord = paste(seed_chr, seed_pos, sep = ":"))

ld_gwas <- distinct(left_join(gwas_efo, ld_snps, by = "seed_snp")) # duplicates likely due to allele switches
# 739 SNPs not run in LD gwas (these were likely ones which had no 1000G match or were monallelic in the chosen population, as per LDlink error message)
check <- anti_join(gwas_efo, ld_snps, by = "seed_snp") ; length(unique(check$seed_snp))

# So bc some snps don't have an LD coord (the ones in check), but we match on that later,
# We should fill in LD_coord with seed_snp coord for missing snps instead!

head(ld_gwas[is.na(rsid_ld),])

# Fill in missing data
ld_gwas[is.na(rsid_ld), rsid_ld := seed_snp]
ld_gwas[is.na(coord_ld), coord_ld := seed_coord] 

# So Far:
length(unique(ld_gwas$seed_snp)); length(unique(ld_gwas$rsid_ld))
# 2901; # 29,396 LD snps


### Join eQTLs to LD_GWAS table ------------------------------------------------

# Join by coord :) # only keep gwas entries which have a SNP in eQTL eGenes table and vice versa
ld_gwas <- rename(ld_gwas, coord = coord_ld) # note LD coords include a seed snp entry
eGenes$coord <- paste(eGenes$chr, eGenes$pos, sep = ":")
ld_gwas_qtl <- clean_cols(inner_join(eGenes, ld_gwas, by = "coord")) %>% distinct()

length(unique(ld_gwas_qtl$disease_trait))
# 376 GWAS eGenes across 4691 loci, 104 disease traits, 37 efo

# Summary

table(ld_gwas_qtl$cell_type)
t_gwas <- filter(ld_gwas_qtl, cell_type == "CD8_T_cells" | cell_type == "CD4_T_cells") 
b_gwas <- filter(ld_gwas_qtl, cell_type == "B_cells") 
m_gwas <- filter(ld_gwas_qtl, cell_type == "Monocytes_Macrophages") 
n_gwas <- filter(ld_gwas_qtl, cell_type == "Neutrophils") 


top_gwas <- ld_gwas_qtl %>% group_by(cell_type, snp) %>% slice_max(p_bf, with_ties = FALSE) 
table(top_gwas$cell_type) 

# Tidy colnames
ld_gwas_qtl <- ld_gwas_qtl %>%
  rename(eQTL_snp=snp, eQTL_type = e_qtl, coord_ld= coord, reported_genes = reported_gene_s, context = contt,
         gwas_rsid = seed_snp, gwas_coord = seed_coord, gwas_study = study) %>%
  select(-m_ua, -m_ub, -chr, -pos, -ref, -alt, -mEff, -is_eGene, -date_added_to_catalog, -first_author, -link, -merged,
         -cnv, -mapped_trait, -mapped_trait_uri, -platform_snps_passing_qc)

column_names <- c("eQTL_snp", "cell_type", "gene_name", "eta",
                  "lrt_eqtl", "p_bf",
                  "gwas_study", "disease_trait","initial_sample_size",
                  "gwas_rsid", "gwas_coord", "region",
                  "reported_genes", "mapped_gene", "snp_gene_ids","upstream_gene_id",
                  "downstream_gene_id", "upstream_gene_distance",
                  "downstream_gene_distance", "strongest_snp_risk_allele", 
                  "context", "intergenic", "risk_allele_frequency",
                  "p_value", "pvalue_mlog", "p_value_tt", "or_or_beta",
                  "x95_percent_ci_tt",
                  "efo_term", "parent_term",
                  "replication_sample_size",
                  "genotyping_technology", "study_accession",
                  "pubmedid", "date", "journal", 
                  "rsid_ld", "coord_ld","alleles", "forg_edb",
                  "regulome_db","function")

ld_gwas_qtl <- ld_gwas_qtl %>% select(column_names) %>% rename(p_adj = p_bf)

rm(zero_ps, files, wd, plotdir, base_dir, gene_sum, data_list)
# save.image(file = "ld_gwas_qtl.RData")


# eQTL gencode annotation ------------------------------------------------------

# Convert SNP data to GRanges
snp_ranges <- GRanges(
  seqnames = eGenes$chr,
  ranges = IRanges(start = eGenes$pos, end = eGenes$pos),
  strand = rep("*", nrow(eGenes))
)

# Get gencode annotation
gtf_fn <- paste0(in_dir, "/SNEEP/resources/gencode.v43.annotation.gtf.gz")
txdb <- GenomicFeatures::makeTxDbFromGFF(file = gtf_fn, format = "gtf")

columns(txdb)
seqlevels(txdb)
seqlevels(txdb) <- paste0(rep("chr"), 1:22)
txdb

# Make genomic ranges objects for each feature
exons <- exonsBy(txdb, by = "tx")
introns <- intronsByTranscript(txdb, use.names = TRUE)
utr5 <- fiveUTRsByTranscript(txdb)
utr3 <- threeUTRsByTranscript(txdb)
upstream <- flank(transcripts(txdb), width = 2000, both = FALSE)
downstream <- flank(transcripts(txdb), width = 2000, both = FALSE, start = FALSE)
promoters <- promot

features <- data.frame(
  exons = rep(0, nrow(eGenes)),
  introns = rep(0, nrow(eGenes)),
  utr5 = rep(0, nrow(eGenes)),
  utr3 = rep(0, nrow(eGenes)),
  upstream = rep(0, nrow(eGenes)),
  downstream = rep(0, nrow(eGenes))
)

features$exons <- overlapsAny(snp_ranges, unlist(exons))
features$introns <- overlapsAny(snp_ranges, unlist(introns))
features$utr5 <- overlapsAny(snp_ranges, unlist(utr5))
features$utr3 <- overlapsAny(snp_ranges, unlist(utr3))
features$upstream <- overlapsAny(snp_ranges, upstream)
features$downstream <- overlapsAny(snp_ranges, downstream)

str(features)
table(features$exons)

# add ncRNA
all_transcripts <- transcripts(txdb, columns = c("gene_id", "tx_id", "tx_type"))
# Filter for non-coding RNAs (adjust the types according to your annotations)
nc_transcripts <- subset(all_transcripts, tx_type %in% c("ncRNA", "lncRNA", "miRNA"))
# Get exons and introns for ncRNAs
nc_exons <- exons(txdb, filter = nc_transcripts)
nc_introns <- unlist(intronsByTranscript(txdb,use.names=TRUE))
nc_introns <- nc_introns[nc_introns$tx_id %in% nc_transcripts$tx_id]


features <- features %>%
  mutate(feature_type = case_when(
    utr3       == TRUE ~ "utr3",
    utr5       == TRUE ~ "utr5",
    exons      == TRUE ~ "exons",
    upstream   == TRUE ~ "upstream",
    downstream == TRUE ~ "downstream",
    introns    == TRUE ~ "introns",
    TRUE               ~ "intergenic"
  ))

table(features$feature_type)


eGenes$feature_type <- features$feature_type
top_gwas <- eGenes %>% group_by(snp) %>% slice_max(p_bf, with_ties = FALSE) %>% mutate(study = "All eQTLs")

cvd_genes <- filter(eGenes, snp %in% ld_gwas_qtl$eQTL_snp)
top_cvd <- cvd_genes %>% group_by(snp) %>% slice_max(p_bf, with_ties = FALSE) %>% mutate(study = "CVD eQTL GWAS")

feature_counts <- bind_rows(top_gwas, top_cvd)
feature_counts <- as.data.frame(table(feature_counts$feature_type, feature_counts$study))
names(feature_counts) <- c("feature_type", "study", "counts")

# Add proportion calcs
feature_counts <- feature_counts %>%
  group_by(study) %>%
  mutate(total_counts = sum(counts)) %>%
  ungroup() %>%
  mutate(proportion = counts / total_counts) %>% select(-total_counts)

feature_counts$feature_type <- as.factor(feature_counts$feature_type)
feature_counts$feature_type <- factor(feature_counts$feature_type , levels = c("introns", "intergenic", "upstream", "downstream", "utr3", "utr5", "exons"))

# Plot feature proportions
# Plotting the proportion data
ggplot(feature_counts, aes(x = feature_type, y = proportion, fill = feature_type)) +
  geom_bar(stat = "identity", color = "black", position = "dodge") +
  theme_minimal() +
  labs(x = "Feature Type", y = "Proportion", title = "Proportion of Genomic Features by Study") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~study)


# Plotting feature counts
ggplot(feature_counts, aes(x = feature_type, y = counts, fill = feature_type)) +
  geom_bar(stat = "identity", color = "black") + 
  theme_minimal() +
  labs(x = "Feature", y = "Count", title = "Genomic Feature Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~study)

# Permute to get significance of these findings

# Distance to TSS --------------------------------------------------------------

tss <- transcriptsBy(txdb, by="gene")
tss <- resize(tss, width=1, fix="start")
nearest_tss <- as.data.frame(distanceToNearest(snp_ranges, unlist(tss)))
nearest_tss <- nearest_tss %>% rename(tss_distance = distance) %>% select(tss_distance)

eGenes <- bind_cols(eGenes, nearest_tss)
top_gwas <- eGenes %>% group_by(snp) %>% slice_max(p_bf, with_ties = FALSE) %>% mutate(study = "All eQTLs")

cvd_genes <- filter(eGenes, snp %in% ld_gwas_qtl$eQTL_snp)
top_cvd <- cvd_genes %>% group_by(snp) %>% slice_max(p_bf, with_ties = FALSE) %>% mutate(study = "CVD eQTL GWAS")


# Plot histogram of distances
study_split <- bind_rows(top_gwas, top_cvd)

study_split %>%
  ggplot(aes(x=tss_distance, fill=study)) +
  geom_density(alpha = 0.5) +
  xlim(0, 1e5) +
  theme_bw() + facet_wrap(~study)

# save.image(file = "ld_gwas_qtl.RData")

#### Load data -----------------------------------------------------------------

# Format tables

load(file = "ld_gwas_qtl.RData")

eGenes <- select(eGenes, -genomic_context, -e_qtl, -is_eGene) # All are cis-eQTLs
eGenes <- rename(eGenes, ensembl_gene_id = gene_name)

# Add gene annotation data -----------------------------------------------------

# Add gene hgnc_name
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

genes <- getBM(
  attributes=c('ensembl_gene_id', 'hgnc_symbol','description','gene_biotype'),
  filters=c('ensembl_gene_id'),
  values = list(ensembl_gene_id=unique(eGenes$ensembl_gene_id)),
  mart = ensembl)

genes <- distinct(genes)
genes <- genes %>% mutate(is_dup = stri_duplicated(ensembl_gene_id)) %>% filter(is_dup != TRUE) %>% select(-is_dup)

eqtl_results <- left_join(eGenes, genes, by = "ensembl_gene_id") # 666 no info eGenes
dim(eGenes) ; dim(eqtl_results)

eqtl_results <- eqtl_results %>% select(snp, cell_type, ensembl_gene_id, hgnc_symbol, description, gene_biotype, m_ua, m_ub, eta, lrt_eqtl, p_nom, mEff, p_bf, feature_type, tss_distance) %>%
  rename(snp_genomic_context = feature_type,
         snp_tss_distance = tss_distance,
         p_adj = p_bf,
         eQTL_snp = snp) 

eqtl_results <- eqtl_results %>% mutate(cvd_gwas_qtl=ifelse(eQTL_snp %in% ld_gwas_qtl$eQTL_snp, TRUE, FALSE))
table(eqtl_results$cvd_gwas_qtl)

# save.image(file = "ld_gwas_qtl.RData")

write.csv(eqtl_results, file = "cseqtl_results/summary_tables/eqtl_results.csv", row.names = F)
write.csv(ld_gwas_qtl, file = "cseqtl_results/summary_tables/eqtl_gwas_results.csv", row.names = F)
write.csv(sig_sneep, file = "cseqtl_results/summary_tables/sneep_results.csv", row.names = F) # in sneep_results.R

# Add snp annotation data ------------------------------------------------------

snpmart <- useEnsembl(biomart = "snp", dataset="hsapiens_snp")
                      
# snpinfo <-  getBM(attributes=c('refsnp_id', 'chr_name','chrom_start', 'allele', 'ensembl_gene_name',
#                                   'minor_allele', 'minor_allele_freq',
#                                   'clinical_significance', 'consequence_type_tv',
#                                   'polyphen_prediction', 'polyphen_score',
#                                   'reg_feature_stable_id', 'reg_consequence_types',
#                                   'motif_consequence_types'),
#                   filters = c('ensembl_gene', 'start'),
#                   values = list(ensembl_gene=unique(eGenes$ensembl_gene_id), start=unique(eGenes$pos)),
#                   mart = snpmart)


snps <- unique(eGenes$pos)
chroms <- unique(eGenes$chr) %>% str_remove_all("chr")
chroms <- c("20") %>% as.numeric()
snps <- filter(eGenes, chr == "chr20") %>% select(pos, ensembl_gene_id)
snp_pos <- unique(snps$pos)
snp_genes <- unique(snps$ensembl_gene_id)[1:2]


snpinfo1 <-  getBM(attributes=c('refsnp_id', 'chr_name','chrom_start', 'allele', 'ensembl_gene_name'),
                  filters = c('ensembl_gene'),
                  values = list(ensembl_gene=snp_genes),
                  mart = snpmart)

# Annovar ----------------------------------------------------------------------

# Select relevant columns and rename them appropriately for ANNOVAR
# 1 based
annovar_data <- eGenes %>% mutate(start = pos, end = pos) %>%
  select(chr, start, end, ref, alt)

# Write the data frame to a text file in tab-separated format
write.table(annovar_data, "annovar/annovar_input.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

annovar_gene <- fread("annovar/annotated.hg38_multianno.txt")
annovar_snp <- fread("annovar/annotated_snp.hg38_multianno.txt") 

wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz


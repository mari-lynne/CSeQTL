# Set up -----------------------------------------------------------------------

library(data.table)
library(parallel)

in_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased/hap"
out_dir=file.path(in_dir, "per_gene")

ref_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/resources/eqtl_500kb22.txt"
ref_data <- fread(ref_dir, header = FALSE)
colnames(ref_data ) <- c("gene_id", "chrom", "start", "end")

# Get sample list --------------------------------------------------------------

rna_ids <- read.csv(file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/metadata/merged/cseqtl_ids_rna_order.csv")

# Function ---------------------------------------------------------------------

# Set up parallel processing ---------------------------------------------------
num_cores <- detectCores() - 1  # Use one less core than available

test_ref <- ref_data[1:34, ]

mclapply(1:nrow(ref_data), function(i) {
  process_gene(test_ref[i, ], in_dir, out_dir, id_dir)
}, mc.cores = num_cores)


# Check ------------------------------------------------------------------------

hap_files <- list.files(in_dir, pattern = "_hap.txt", full.names = TRUE)
hap_geno_ids <- gsub("_hap.txt", "", basename(hap_files))
hap_files <- hap_files[hap_geno_ids %in% rna_ids$geno_id]

process_gene <- function(gene_row, hap_files, output_dir) {  # Opening function block
  gene_id <- gene_row$gene_id
  chrom <- gene_row$chrom
  start <- gene_row$start
  end <- gene_row$end
  
  # Initialize an empty data table to store combined data
  combined_data <- data.table()
  
  # Loop through all hap files
  for (hap_file in hap_files) {  # Opening loop block
    hap_data <- fread(hap_file, header = FALSE, fill = TRUE)
    
    # Remove any rows where the first column is NA
    hap_data <- hap_data[!is.na(V1)]  
    
    # Assign column names
    colnames(hap_data) <- c("chrom", "pos", "ref", "alt")
    
    # Subset hap data based on gene coordinates
    subset_data <- hap_data[chrom == as.character(chrom) & pos >= start & pos <= end,]
    
    # Add sample identifier
    sample_id <- gsub("_hap.txt", "", basename(hap_file))
    subset_data[, sample_id := sample_id]
    
    # Bind to combined data table
    combined_data <- rbind(combined_data, subset_data, fill = TRUE)
  }  # Closing loop block
  
  # Write combined data to output file
  output_file <- file.path(output_dir, paste0(gene_id, "_genotype.txt"))
  fwrite(combined_data, output_file, sep = "\t", col.names = TRUE, row.names = FALSE)
}  # Closing function block



test_ref <- ref_data[1:4, ]
small_hap_files <- hap_files[1]  # Use only the first haplotype file
test_result <- process_gene(test_ref[2, ], small_hap_files, out_dir)




# Set up parallel processing ---------------------------------------------------
num_cores <- detectCores() - 1  # Use one less core than available

test_ref <- ref_data[1:4, ]

mclapply(1:nrow(ref_data), function(i) {
  process_gene(test_ref[i, ], hap_files, out_dir)
}, mc.cores = num_cores)

test_result <- lapply(1:nrow(test_ref), function(i) {
  process_gene(test_ref[i, ], hap_files, out_dir)
})



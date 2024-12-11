# Aims

# Filter TReC and ASReC files to exclude transcripts with low counts as identified in rna_preprocess.R
# Run on FhR with 8 CPUs

library(data.table)
library(stringr)

# Set the path to the input directory
input_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/rnaseq/ASE"
setwd(input_dir)

# Set the path to the output directory
output_dir <- file.path(input_dir, "qc_genes")

# Set the path to the gene filter file
gene_filter_file <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/ASE/redo/rna_qc_filter_genes.txt"

# Create the output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE)

# Read the gene filter file into a vector
filter_genes <- readLines(gene_filter_file)

# Loop -------------------------------------------------------------------------

start_time <- Sys.time()
# Loop through each file in the input directory
for (file in list.files(path = input_dir, pattern = "-output.trecase.txt", full.names = TRUE)) {
  
  # Extract sample_id from the file name
  sample_id <- sub("-output.trecase.txt", "", basename(file))
  
  setDTthreads(0L)
  # Read the input file using data.table
  data <- fread(file) %>% setnames("V1", "ensemble_gene_id")
  data$ensemble_gene_id <- str_remove(data$ensemble_gene_id, "\\..*")
  
  # Set the ensemble_gene_id as the key for faster filtering
  setkey(data, ensemble_gene_id)
  
  # Filter rows based on gene ID
  filtered_data <- data[filter_genes, nomatch = 0]
  # Remove dup ids (~20 genes)
  filtered_data <- filtered_data[!duplicated(filtered_data$ensemble_gene_id), ]
  
  # Write the filtered data to the output file
  write.table(filtered_data, file.path(output_dir, paste0(sample_id, "-output.trecase.txt")), sep = "\t", row.names = FALSE, quote = FALSE)
}

end_time <- Sys.time()
time_taken <- format(end_time - start_time, format = "%H:%M:%S")
print(paste("Time taken to filter genes:", time_taken))

# Empty files ------------------------------------------------------------------

# Make empty text files to fill:
output_dir2 <- file.path(output_dir, "/per_gene")
  
# Create empty text files for each filter gene
for (gene_id in filter_genes[-1]) {
  # Define the filename
  filename <- paste0(gene_id, "_hap_counts.txt")
  
  # Create an empty file
  file.create(file.path(output_dir2, filename))
}



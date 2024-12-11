#!/usr/bin/Rscript

library(data.table)
library(stringr)
library(parallel) 

# Aims
# - Aggregate all aseq files
# - Split and save per gene files

# Set up =======================================================================

# INPUT
# Set the path to the directory containing input sample files
input_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/rnaseq/ASE/qc_genes"

# Get a list of all sample files in the directory
aseq_files <- list.files(path = input_dir, pattern = "-output.trecase.txt", full.names = TRUE)

# Get the number of cores available from SLURM allocation
num_cores <- 36  # As specified in the SLURM script

# OUTPUT
output_dir <- file.path(input_dir, "per_gene")

if (dir.exists(output_dir) == FALSE) {
  dir.create(output_dir, showWarnings = FALSE)
}

# Function to process each aseq file (per sample)
process_aseq_file <- function(aseq_file) {
  # Extract sample_id from the file name
  sample_id <- sub("-output.trecase.txt", "", basename(aseq_file))
  # Read the aseq file using data.table
  aseq_data <- fread(aseq_file, skip = 1)  # Skip the header
  # Add a column for sample_id
  aseq_data$sample_id <- sample_id
  return(aseq_data)
}

# Aggregate all trecase files using parallel processing =========================

setDTthreads(0L) 
agg_gene_counts <- data.table()

start_time <- Sys.time()
# Use mclapply for parallel processing
aseq_data_list <- mclapply(aseq_files, process_aseq_file, mc.cores = num_cores)
# Combine all the data tables into one
agg_gene_counts <- rbindlist(aseq_data_list, use.names = TRUE, fill = TRUE)
end_time <- Sys.time()

time_taken <- format(end_time - start_time, format = "%H:%M:%S")
print(paste("Time taken to aggregate aseq per gene files:", time_taken))

# Split and save per gene files ================================================

setDTthreads(0L)
transcripts <- unique(agg_gene_counts$V1)
start_time_split <- Sys.time()

# Use mclapply for parallelization
output_files <- mclapply(transcripts, function(transcript) {
  # Subset data for the current transcript
  transcript_data <- agg_gene_counts[V1 == transcript]
  
  # tidy df/names
  transcript_data <- transcript_data[, .(sample_id = sample_id, trec = V2, hap1 = V3, hap2 = V4, hapN = V5)]
  
  # Set the path for the output gene file for the current transcript
  output_file_name <- file.path(output_dir, paste0(transcript, "_hap_counts.txt"))
  
  # Write the data to the output file
  write.table(transcript_data, output_file_name, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Return the output file name for tracking
  return(output_file_name)
}, mc.cores = num_cores)

end_time_split <- Sys.time()
time_taken_split <- format(end_time_split - start_time_split, format = "%H:%M:%S")
print(paste("Time taken to split and save per transcript files:", time_taken_split))

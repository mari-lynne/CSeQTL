# Aims

# Get roadmap chipseq data for relevant cell types
# Use liftover for converting bed files to Hg38
# Test for enrichment of cs-eqtls in cs-eqtl chipseq data

# Roadmap download -------------------------------------------------------------

# Load necessary libraries
library(rvest)
library(R.utils)  # For gunzip function
options(timeout = 400)  # Increase allowed download time

# Define function to download and extract Roadmap files
download_roadmap_files <- function(cell_type_code = "E029", 
                                   out_folder = getwd(), 
                                   url = "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/") {
  
  # Read the HTML content from the webpage
  webpage <- read_html(url)
  
  # Extract the links to the files
  file_links <- webpage %>%
    html_nodes("a") %>%
    html_attr("href")
  
  # Filter for files starting with the specified cell type code
  selected_files <- file_links[grepl(paste0("^", cell_type_code), file_links)]
  
  # Full URLs for the files
  file_urls <- paste0(url, selected_files)
  
  # Create output folder if it does not exist
  if (!dir.exists(out_folder)) {
    dir.create(out_folder, recursive = TRUE)
  }
  
  # Download and unzip each file
  for (file_url in file_urls) {
    # Extract the file name
    file_name <- basename(file_url)
    
    # Define full path for the destination file
    dest_file <- file.path(out_folder, file_name)
    
    # Download the file
    download.file(file_url, destfile = dest_file, mode = "wb")
    cat("Downloaded:", file_name, "\n")
    
    # Unzip the file if it ends with .gz
    if (grepl("\\.gz$", file_name)) {
      # Unzip the file
      unzipped_file <- sub("\\.gz$", "", dest_file)  # Remove .gz extension
      gunzip(dest_file, destname = unzipped_file, overwrite = TRUE)
      cat("Unzipped:", unzipped_file, "\n")
    }
  }
}


# Example Usage
# download_roadmap_files(
#   cell_type_code = "E029",
#   out_folder = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/resources/chipseq/monocytes",
#   url = "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/"
# )
# In future add option for selecting narrow or broad URL


download_roadmap_files(
  cell_type_code = "E034",
  out_folder = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/resources/chipseq/t_cells",
  url = "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/"
)


# Roadmap liftover -------------------------------------------------------------

# Update peak bed files to Hg38 Format

# /fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/SCT/eqtl/Downstream
# ./liftover.sh

# Prep SNP BED files -----------------------------------------------------------

# Update SNP results files to BED format for BEDTools intersect

load(file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT/eqtl_multiple_testing.RData")

setwd("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/snp_nexus/output")

options(scipen = 999)

# Background SNP file:
write.table(file = "background_snps_38.bed",
  data.frame(chr = eqtls$chr, start = eqtls$pos - 1, end = eqtls$pos, snp = eqtls$snp), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

# Save CS- SNP files:

# Loop thru list of data frames (e.g., B_cells, T_cells, etc.)
lapply(names(results), function(name) {
  df <- results[[name]]
  
  # Create and save BED format in one step
  write.table(file = paste0(name, "_snps.bed"),
    data.frame(chr = df$chr, start = df$pos - 1, end = df$pos, snp = df$snp),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
  )
})

# In gui folders move files to subdirs


# Find SNP/Epigenetic Feature Overlaps  ----------------------------------------

# Run BED tools intersect
# Will run permutation script in parallel across features per cell type
# Therefore make an array file of features per cell type

# #!/bin/bash
# 
# # Set base directory
# base_dir="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl"
# 
# # Define cell types
# cell_types=("monocytes" "b_cells" "t_cells" "neuts")
# 
# # Loop through each cell type to generate array files
# for cell_type in "${cell_types[@]}"; do
# 
# # Set the peak directory for the current cell type
# peak_dir="${base_dir}/resources/chipseq/${cell_type}/lifted"
# 
# # List all .narrowPeak files and save them to the array file
# array_file="${base_dir}/scripts/SCT/eqtl/Downstream/${cell_type}_array.txt"
# ls "${peak_dir}"/*.narrowPeak > "$array_file"
# 
# echo "Array file created for $cell_type: $array_file"
# 
# done


#  sbatch enrich_permute2.sh b_cell_array.txt

# In terminal
# ml BEDTools
# bedtools intersect -a hglft_genome_219_ebf280.bed -b E032-DNase.hotspot.broad.bed -wa -wb > DNase_overlaps.txt
# bedtools intersect -a hglft_genome_219_ebf280.bed -b E032-DNase.hotspot.broad.bed -wa -wb > DNase_overlaps.txt




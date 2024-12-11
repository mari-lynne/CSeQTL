# data <- fread("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/rnaseq/ASE/NWD106483-output.trecase.txt") %>% setnames("V1", "ensemble_gene_id")
# data$ensemble_gene_id <- str_remove(data$ensemble_gene_id, "\\..*")
# setkey(data, ensemble_gene_id)
# filtered_data <- data[filter_genes, nomatch = 0]
# dups <- filtered_data[duplicated(filtered_data$ensemble_gene_id), ]

library(ggplot2)
library(data.table)
library(stringr)
library(dplyr)
library(tidylog)
library(ggpubr)

# LLS V2 -----------------------------------------------------------------------

sample_id <- "NWD106483"

data2 <- fread("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/rnaseq/ASE/qc_genes/NWD106483-output.trecase.txt") 

data2$ensemble_gene_id <- str_remove(data2$ensemble_gene_id, "\\..*")

# Reshape the data from wide to long format, excluding the ensemble_gene_id
summary_data2  <- data2 %>%
  pivot_longer(cols = -ensemble_gene_id, names_to = "sample", values_to = "count")

# Summarize the total values for each column
summary_data2 <- summary_data2 %>%
  group_by(sample) %>%
  summarize(total = sum(count, na.rm = TRUE))

# Create the bar chart using ggplot2
b <- ggplot(summary_data2, aes(x = sample, y = total)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Total and Haplotype Specific RNA-Seq Counts (LLS V2)",
       x = "Sample",
       y = "Total Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# LLS V1 -----------------------------------------------------------------------

sample_id <- "NWD106483"

data <- fread("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/ASE/redo/NWD106483-output.trecase.txt") %>% setnames("V1", "ensemble_gene_id")
data$ensemble_gene_id <- str_remove(data$ensemble_gene_id, "\\..*")

# Filter for genes in data
data <- data[data$ensemble_gene_id %in% data2$ensemble_gene_id, ]

# Reshape the data from wide to long format, excluding the ensemble_gene_id
summary_data  <- data %>%
  pivot_longer(cols = -ensemble_gene_id, names_to = "sample", values_to = "count")

# Summarize the total values for each column
summary_data <- summary_data %>%
  group_by(sample) %>%
  summarize(total = sum(count, na.rm = TRUE))

# Create the bar chart using ggplot2
a <- ggplot(summary_data, aes(x = sample, y = total)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Total and Haplotype Specific RNA-Seq Counts (LLS V1)",
       x = "Sample",
       y = "Total Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(patchwork)
(a | b)


# Loop all samples (v2) -------------------------------------------------------------

# Step 1: List all files in the directory with the pattern "-output.trecase.txt"
file_list <- list.files(path = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/rnaseq/ASE/qc_genes",
                        pattern = "-output.trecase.txt", full.names = TRUE)

files_ids <- str_extract(basename(file_list), "^[^-]+")
files_ids <- as.numeric(str_extract(files_ids, "\\d+"))

# Step 2: Order the file_list numerically based on the sample IDs
ordered_file_list <- file_list[order(files_ids)]

file_list <- file_list[1:120]

# Step 2: Initialize an empty list to store data
all_data <- list()

# Step 3: Loop through each file, read data, extract sample ID, and store in the list
for (file in file_list) {
  # Extract the sample ID from the filename
  sample_id <- basename(file) %>% str_extract("^[^-]+")
  
  # Read the data
   data <- fread(file) # %>% setnames("V1", "ensemble_gene_id")
  # data$ensemble_gene_id <- str_remove(data$ensemble_gene_id, "\\..*")
  
  # Add the sample ID as a column
  data$sample_id <- sample_id
  
  colnames(data) <- c("ensemble_gene_id", "TReC", "Hap1", "Hap2", "HapN", "sample_id")
  
  # Append the data to the list # index by sample id
  all_data[[sample_id]] <- data
}

# Step 4: Combine all the data into a single data frame
combined_data <- bind_rows(all_data) # 828 files read (i got bored waiting)

# Step 5: Reshape the data from wide to long format
long_data <- combined_data %>%
  pivot_longer(cols = -c(ensemble_gene_id, sample_id), 
               names_to = "category", 
               values_to = "count")

# Step 6: Summarize the total counts for each sample and each category
summary_data <- long_data %>%
  group_by(sample_id, category) %>%
  summarize(total = sum(count, na.rm = TRUE), .groups = 'drop')

# Step 7: Create the plot
a2 <- ggplot(summary_data, aes(x = category, y = total, color = sample_id)) +
  geom_point() +
  theme_minimal() +
  labs(title = "ASeq Mapping LLS (V2)",
       x = "Category",
       y = "Total Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") 

a2

a22 <- summary_data %>% 
  filter(category != "TReC") %>%
  ggboxplot(x = "category", y = "total", 
            color = "category", palette = "jco",
            add = "jitter", 
            ggtheme = theme_minimal()) +
  labs(title = "ASeq Mapping LLS (V2)",
       x = "Category",
       y = "Total Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
a22


filter_genes <- unique(data$ensemble_gene_id)

# Loop all samples (v1) -------------------------------------------------------------

# Step 1: Extract sample IDs from the filenames in the file_list
file_list <- list.files(path = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/ASE/redo/",
                        pattern = "-output.trecase.txt", full.names = TRUE)

files_ids <- str_extract(basename(file_list), "^[^-]+")
files_ids <- as.numeric(str_extract(files_ids, "\\d+"))

# Step 2: Order the file_list numerically based on the sample IDs
ordered_file_list <- file_list[order(files_ids)]

file_list <- ordered_file_list[1:120]

# Instialize empty list
all_data1 <- list()

# Step 3: Loop through each file, read data, extract sample ID, and store in the list
for (file in file_list) {
  # Extract the sample ID from the filename
  sample_id <- basename(file) %>% str_extract("^[^-]+")
  
  # Read the data
  data1 <- fread(file)
  colnames(data1) <- c("ensemble_gene_id", "TReC", "Hap1", "Hap2", "HapN")
  data1$ensemble_gene_id <- str_remove(data1$ensemble_gene_id, "\\..*")
  
  # Add the sample ID as a column
  data1$sample_id <- sample_id
  
  # Filter for genes in filter_genes
  data1 <- data1[data1$ensemble_gene_id %in% filter_genes, ]
  
  # Append the data to the list # index by sample id
  all_data1[[sample_id]] <- data1
}

# Step 4: Combine all the data into a single data frame
combined_data1 <- bind_rows(all_data1) # 828 files read (i got bored waiting)

# Step 5: Reshape the data from wide to long format
long_data1 <- combined_data1 %>%
  pivot_longer(cols = -c(ensemble_gene_id, sample_id), 
               names_to = "category", 
               values_to = "count")

# Step 6: Summarize the total counts for each sample and each category
summary_data1 <- long_data1 %>%
  group_by(sample_id, category) %>%
  summarize(total = sum(count, na.rm = TRUE), .groups = 'drop')  
  
### Plot =----------------------------------------------------------------------

# Step 7: Create the plot
a1 <- ggplot(summary_data1, aes(x = category, y = total, color = sample_id)) +
  geom_point() +
  theme_minimal() +
  labs(title = "ASeq Mapping LLS (V1)",
       x = "Category",
       y = "Total Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") 

a1

a11 <- summary_data %>% 
  filter(category != "TReC") %>%
  ggboxplot(x = "category", y = "total", 
            color = "category", palette = "jco",
            add = "jitter", 
            ggtheme = theme_minimal()) +
  labs(title = "ASeq Mapping LLS (V1)",
       x = "Category",
       y = "Total Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
a11

(a11 | a22)

ggsave(filename = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/plots/merged/aseq_hap_box_lls_comparison.png")

summary_data %>% 
  filter(category != "TReC") %>%
  ggboxplot(x = "category", y = "total", 
            color = "category", palette = "jco",
            add = "jitter", 
            ggtheme = theme_minimal()) +
  labs(title = "ASeq Mapping",
       x = "Category",
       y = "Total Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(filename = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/plots/merged/aseq_hap_boxplot.png")

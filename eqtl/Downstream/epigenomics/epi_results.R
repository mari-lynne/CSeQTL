# Load packages
library(data.table)
library(dplyr)
library(stringr)
library(janitor)
library(tidylog)
library(ggplot2)

# Results ----------------------------------------------------------------------

# CSeQTL Resuls (from mt_nov script.R)
eqtls <- read.csv(file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT/mt_nov_results.csv")

cell_type="neuts"
base_dir="/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/snp_nexus/output"
in_dir=file.path(base_dir, paste0(cell_type, "/peak_overlaps/cseqtl_overlap_snps"))
setwd(in_dir)

plot_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/plots"

# results.txt files have overlapping snps (redo to add exact epi feature TODO)
# results.csv files have permutation results

file_names <- list.files(path=in_dir, all.files = TRUE, full.names = FALSE, no..=TRUE) # Get file names
header <- c("chr", "start", "end", "snp")# add column names

# Read all files, assign column names, and add the "feature" column
myfiles <- lapply(file_names, function(file) {
  # Extract file name without path and extension
  base_name <- tools::file_path_sans_ext(basename(file))
  feature_name <-  str_extract(base_name, "E[0-9]+-[A-Za-z0-9]+")
  
  # Read the file and assign column names
  df <- fread(file, col.names = header)
  
  # Add the "feature" column with the file base name
  df$feature <- feature_name
  df$cell_type <- cell_type
  
  # Return the modified data frame
  return(df)
})

names(myfiles) <- str_extract(file_names, "E[0-9]+-[A-Za-z0-9]+")
epi_features <- rbindlist(myfiles) # list2env((myfiles), envir = .GlobalEnv)

b_epi <- epi_features
t_epi <- epi_features
m_epi <- epi_features
n_epi <- epi_features

epi_features <- bind_rows(b_epi, t_epi, m_epi, n_epi)
epi_features$feature <- str_remove_all(epi_features$feature, "E[0-9]+-")

# Quick summary
table(epi_features$cell_type, epi_features$feature)

# Plot epi feautre data --------------------------------------------------------

epi_features %>%
  filter(feature != "DNase" | feature != "H3K27ac") %>%
  ggplot(aes(x = cell_type, fill = feature))  +
  geom_bar() + 
  labs(x = "Cell Type", y = "Count", fill = "Epigenetic Feature") +
  theme_minimal() +  # A clean minimalistic theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# DNAse- open chromatin
# H3K27ac - enhancers
# H3K4me1- primed enhancers
# H3K4me3-promoters
# H3K36me3-gene bodies
# H3K27me3-polycomb repression
# H3K9me3-heterochromatin

epi_features$transcription_association <- 
  ifelse(epi_features$feature %in% c("DNase", "H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3"),
                                                 "activation", "repression")

# Sort levels:

# Convert 'feature' to a factor with specified levels in the desired order
epi_features$feature <- factor(epi_features$feature,
                               levels = c("DNase", "H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3"))

# Convert 'feature' to a factor with specified levels in the desired order
epi_features$cell_type <- factor(epi_features$cell_type,
                               levels = c("monocytes", "b_cells", "t_cells", "neuts"))

# Filter features not found in all subsets
epi_features2 <- epi_features %>% 
  filter(feature != "DNase" & feature != "H3K27ac") %>% droplevels()

# Calculate proportions
epi_features_prop <- epi_features2 %>%
  group_by(cell_type, feature) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(total = sum(count), proportion = count / total)

# Create a bar chart with proportions
ggplot(epi_features_prop, aes(x = cell_type, y = proportion, fill = feature)) +
  geom_bar(stat = "identity", position = "fill") +  # Use fill to stack proportions to sum to 1
  labs(x = "Cell Type", y = "Proportion", fill = "Feature", title = "Epigenetic Features across Cell Types") + scale_fill_manual(values=c("#FF595E","#FFCA3A","#8AC926","#1982C4","#6A4C93")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", colour = NA)) 

ggsave(file = file.path(plot_dir, "epi_feature_sum.png"))

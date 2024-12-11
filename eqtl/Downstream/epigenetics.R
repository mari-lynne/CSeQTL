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


# Load data --------------------------------------------------------------------

cell_type <- "b_cell"

workdir <- file.path("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/snp_nexus/output", cell_type)
setwd(workdir)

plotdir <- c("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/plots")

# Load eQTL data
load(file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT/eqtl_multiple_testing.RData")

# Load epigenetic data
roadmap <- clean_names(fread("roadmap.txt"))

# Summarize Road Map data ------------------------------------------------------

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


# Cell type checks -------------------------------------------------------------

b_epi <- filter(roadmap, str_detect(epigenome, "B cell"))
# 1630 features in B-cells matching GWAS SNPs

# To permute I would need to upload random (both non eQTL and those from CSeQTL data) SNPs to SNP nexus each time and see how many B-cell SNPs I get

simple <- b_epi %>% select(variation_id, chromosome, feature_type)
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
  labs(fill = "Feature Type", title = "Distribution of B-cell Histone Markers")




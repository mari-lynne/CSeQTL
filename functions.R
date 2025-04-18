# Functions 

# TODO

# clean_left_join function that drops dup .y vars and cleans colnames by removing ".y"
# wrapper function to pate dirs and file names when write.csv or write.table

# Helper functions -------------------------------------------------------------

# use with data %>% select(where(not_na))
not_all_na <- function(x) all(!is.na(x))
not_any_na <- function(x) any(!is.na(x)) 

# Not in vector
'%!in%' <- function(x,y)!('%in%'(x,y))

# TPM - transcripts per million

tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

clean_cols <- function(df) {
  # Tidy duplicate columns
  colnames(df) <- str_remove(colnames(df), ".x")
  # Remove duplicate column
  df <- df[ ,str_detect(colnames(df), ".y")==FALSE]
  # Remove columns with all NA values
  df <- df[, colSums(is.na(df)) != nrow(df)]
  return(df)
}

# Contrast matrix function -----------------------------------------------------

#/ Define groups and baseline to make contrasts
make_contrasts <- function (group,
                            control,
                            delim = "_vs_",
                            des_mat) {
  
  suppressMessages(require(limma))
  
  # Checks
  if (is.null(group))
    stop("Error: group arg is missing")
  
  #/ Ensure unique group levels
  group <- sort(unique(as.character(group)))
  
  #/ Write code for limma by pasting groups
  # IF control var is present, compare all groups to control
  # ELSE make all comparisons using combn function
  
  if (!missing(control)) {
    combo <- paste0(group, "-",  control)
  } else{
    combo <- combn(group,2,
                   FUN = function(x) {
                     paste0(x[1], "-", x[2])
                   }
    )
  }
  
  #/ Make contrasts
  #/ Option to provide design matrix, or use group levels
  if (!missing(des_mat)) {
    contrasts <-
      limma::makeContrasts(contrasts = combo, levels = colnames(des_mat))
  } else{
    contrasts <- limma::makeContrasts(contrasts = combo, levels = group)
    message("No Design Matrix provided, using only defined contrasts for matrix")
  }
  colnames(contrasts) <- gsub("-", delim, colnames(contrasts))
  return(contrasts)
}

# Todo:
# Levels need to be design matrix, poss modify to have group or design option later
# Only do lmfit step if design and fit are supplied

# Contrasts and ebyaes function ------------------------------------------------

contrast_2_lm <-
  function(group,
           control,
           delim = "_vs_",
           des_mat,
           efit,
           topTab = "TRUE") {
    # Define groups and baseline to make contrasts
    # Input design and previous model fit model for it to work
    # Toptab gives option to output results
    
    #/ Checks
    suppressMessages(require(limma))
    if (is.null(group))
      stop("Error: group arg is missing")
    
    #/ Ensure unique group levels
    group <- sort(unique(as.character(group)))
    
    #/ Define contrasts
    if (!missing(control)) {
      # compare to control
      combo <- paste0(group, "-",  control)
    } else{
      # make all pairwise comparisons
      combo <-
        combn(group,2,
              FUN = function(x) {
                paste0(x[1], "-", x[2])
              }
        )
    }
    
    #/ Make contrast matrix
    if (!missing(des_mat)) {
      contrasts <-
        limma::makeContrasts(contrasts = combo, levels = colnames(design))
    } else{
      #No design + no efit
      contrasts <- limma::makeContrasts(contrasts = combo, levels = group)
      message("No Design Matrix provided, using only defined contrasts for matrix")
    }
    colnames(contrasts) <- gsub("-", delim, colnames(contrasts))
    
    
    #/ Model fit and DEG results
    
    # Requires efit and dmat args to be supplied, else if error message
    # Option to
    if (!missing(efit) &
        !missing(des_mat) & topTab == "TRUE") {
      fit2 <- contrasts.fit(efit, contrasts)
      fit2 <- eBayes(fit2)
      message("Performing ebayes fit of linear model")
      top_results <- list() #TopTable Results
      for (i in colnames(contrasts)) {
        top_results[[i]] <- topTable(fit2, coef = i, number = Inf)
        limma_list <-
          list(contrasts = contrasts,
               fit2 = fit2,
               top_results = top_results)
      }
      message(
        "Contrast matrix (contrasts), ebayes (fit2), top DEGs (top_results) saved in list\n
          Subset list with either $ or [[]] for results"
      )
    } else if (!missing(efit) &
               !missing(des_mat) & topTab != "TRUE") {
      # Without Top table option
      fit2 <- contrasts.fit(efit, contrasts)
      fit2 <- eBayes(fit2)
      message("Performing ebayes fit of linear model")
      limma_list <- list(contrasts = contrasts, fit2 = fit2)
    } else if (missing(efit) | missing(des_mat)) {
      # No top tab or design matrix supplied: ERROR
      limma_list <- list(contrasts)
      warning("No linear model or design matrix supplied, returning contrast matrix only")
    }
    return(limma_list)
  }

# Volcano plot function --------------------------------------------------------

plot_vol <- function(deg,
                     p = 0.05,
                     FC = 0,
                     lab_type = "top",
                     genes,
                     top = 20,
                     title = "",
                     alpha = 0.98,
                     colours = c("#800000", "#a50000", "#ef5a3a", "orange", "yellow")) {
  #/ 1) Checks
  if (missing(deg)) {
    stop("Error: deg arg is missing. Please provide a toptable data frame")
  }
  if (!missing(lab_type) & lab_type != c('top') & missing(genes)) {
    stop("Error: label requires character vector with selected genes")
  }
  
  #/ 2) Define sig values and thresholds
  log_p <- -log10(p)
  # Adj p values #need to convert adjusted p-val to un-adjusted for plot y-scale
  values <- seq(0.050, 0.051, by = 0.00001)
  deg$adj.P.Val <-
    round(deg$adj.P.Val, 5) # This is such a janky way I apologize
  adj.p <- deg[(deg$adj.P.Val %in% (values)), c("P.Value")]
  adj.p <-  min(adj.p)
  log_adj <- -log10(adj.p)
  
  #/ 3) Set up DEG table
  deg <- deg %>%
    mutate(
      reg =
        case_when(
          deg$logFC >= FC & deg$adj.P.Val <= p ~ "Sig Adj. P <0.05",
          deg$logFC <= FC &
            deg$adj.P.Val <= p ~ "Sig Adj. P <0.05",
          deg$logFC >= FC & deg$P.Value <= p ~ "Sig P <0.05",
          deg$logFC <= FC & deg$P.Value <= p ~ "Sig P <0.05",
          abs(deg$logFC) <= FC &
            deg$adj.P.Val >= p ~ "No Change",
          abs(deg$logFC) <= FC &
            deg$adj.P.Val <= p ~ "No Change",
          abs(deg$logFC) > FC &
            deg$adj.P.Val > p ~ "No Change"
        )
    ) %>%
    mutate(reg =
             factor(reg, levels =
                      c(
                        "Sig Adj. P <0.05", "Sig P <0.05", "No Change"
                      )))
  
  #/ 4) Define labels
  
  if (is.null(lab_type)) {
    #No entry for lab list, currently there's no default arg
    gene_label <- c("")
    lab_data <-
      NULL #Empty dataframe, need to work out how not to error
    warning("No genes highlighted")
  } else if (lab_type == "adj.sig") {
    gene_label <- genes
    lab_data <-
      deg[(deg$reg == "Sig Adj. P <0.05") & (deg$gene_name %in% genes), ]
  } else if (lab_type == "sig") {
    gene_label <- genes
    lab_data <-
      deg[(deg$reg == "Sig P <0.05"| deg$reg =="Sig Adj. P <0.05") & (deg$gene_name %in% genes), ]
  } else if (lab_type == "ns") {
    gene_label <- genes
    lab_data <- deg[(deg$gene_name %in% genes), ]
  } else if (lab_type == "top") {
    lab_data <- slice_min(deg, adj.P.Val, n = top)
    gene_label <- lab_data$gene_name
  }
  
  # Default option, then plot for if FC = "shade"
  
  # 5) Plot Volcano
  
  vol <-
    deg %>% ggplot(aes(
      x = logFC,
      y = -log10(P.Value),
      label = gene_name
    )) +
    geom_point(aes(color = P.Value, alpha = alpha)) +
    labs(title = title) +
    theme_minimal() + theme(legend.position = "none") +
    geom_hline(yintercept = log_p,
               linetype = 2.5,
               alpha = 0.7) +
    geom_hline(yintercept = log_adj,
               linetype = 2.5,
               alpha = 0.7) +
    geom_label_repel(
      data = lab_data,
      size = 3.5,
      direction = "both",
      nudge_y = 1.6,
      nudge_x = 0.1,
      angle = 70,
      vjust = 0,
      segment.size = 0.5,
      segment.color = "#331002",
      fill = "#f7f7f5"
    ) +
    scale_color_gradientn(colours = colours,
                          values = c(0, adj.p, p, 1))
  return(vol)
} 


# Read geno function -----------------------------------------------------------

read_omic <- function(name = "", dir = "", ...) {
  # Load required packages
  if (!requireNamespace("stringr", quietly = TRUE)) {
    install.packages("stringr")
  }
  library(stringr)
  
  # Use the current working directory if none is specified
  if (is.null(dir) || dir == "") {
    dir <- getwd()
  }
  
  # Check file type and read accordingly
  if (str_ends(name, "\\.csv")) {
    data <- read.csv(file.path(dir, name))
    
  } else if (str_detect(name, "\\.gct")) {
    data <- read.delim(file.path(dir, name), skip = 2)
    
  } else if (str_ends(name, "\\.rda") || str_ends(name, "\\.RData") || str_ends(name, "\\.Rdata")) {
    data <- load(file.path(dir, name))
    
  } else if (str_ends(name, "\\.rds")) {
    data <- readRDS(file.path(dir, name))    
    
  } else if (str_ends(name, "\\.ods")) {
    if (!requireNamespace("readODS", quietly = TRUE)) {
      install.packages("readODS")
    }
    library(readODS)
    data <- read_ods(file.path(dir, name))
    
  } else {
    stop("Please specify a valid data type: .csv, .gct, .rda, .rds, .RData, .ods")
  }
  
  return(data)
}

# Plot PCAs --------------------------------------------------------------------

# Function to create a 3 x 3 grid of adjacent PC scatter plots
plot_pca_grid <- function(data, num_components = 10, var = "", pal = NULL, title = "",
                          save = FALSE, dir = getwd(), file_name = "plot",
                          width = 11, height = 8) {
  plots <- list()
  
  for (i in 1:(num_components - 1)) {
    pc1_index <- i
    pc2_index <- i + 1
    
    plot_data <- data[, c(paste0("PC", pc1_index), paste0("PC", pc2_index), var)]
    
    plot <- ggplot(plot_data,
                   aes_string(
                     x = paste0("PC", pc1_index),
                     y = paste0("PC", pc2_index),
                     color = var
                   )) +
      geom_point(alpha = 0.8, size = 1.2) +
      theme_bw()
    
    # Make palette if not supplied and using factor var
    if (is.null(pal) && !is.null(var)) {
      if (is.factor(data[[var]])) {
        new_pal <- viridis(n = length(levels(data[[var]])))
        plot <- plot + scale_color_manual(values = new_pal)
      } else { # If numeric just use viridis scale
        plot <- plot + scale_color_viridis()
      }
    } else if (!is.null(pal) && !is.null(var)){ # If palette supplied 
      plot <- plot + scale_color_manual(values = pal)
    }
    
    plots[[i]] <- plot 
  }
  
  grid_matrix <-
    wrap_plots(plots, guides = "collect") + plot_annotation(title = title)
  
  if (save == TRUE) {
    ggsave(filename = paste0(dir, "/", file_name, ".pdf"), plot = grid_matrix, device = "pdf", width = width, height = height)
    ggsave(filename = paste0(dir, "/", file_name, ".png"), plot = grid_matrix, device = "png", width = width, height = height)
  }
  return(grid_matrix)
}

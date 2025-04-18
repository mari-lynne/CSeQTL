library(ggplot2)
library(dplyr)

# Function to create a 3 x 3 grid of adjacent PC scatter plots
# TODO add colour pal functions 
create_pca_grid <- function(data, num_components = 10, var = "ethnic") {
  plots <- list()
  
  for (i in 1:(num_components - 1)) {
    pc1_index <- i
    pc2_index <- i + 1
    
    plot_data <- data %>% select(.data[[paste0("PC", pc1_index)]], .data[[paste0("PC", pc2_index)]], var)
    
    plot <- plot_data %>%
      ggplot(aes(.data[[paste0("PC", pc1_index)]], .data[[paste0("PC", pc2_index)]], color = ethnic)) +
      geom_point(alpha = 0.8, size = 1.2) +
      scale_color_manual(values = pal) +
      theme_bw() +
      labs(color = "Ethnicity") +
      ggtitle(paste0("PC", pc1_index, " vs PC", pc2_index))
    
    plots[[i]] <- plot
  }
  
  # Combine the plots into a 3 x 3 grid
  grid_matrix <- wrap_plots(plots, guides = "collect")
  
  return(grid_matrix)
}

# Usage
# pca_grid_plot <- create_pca_grid(pca_df)
# print(pca_grid_plot)

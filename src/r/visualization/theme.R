#' iHMP Analysis Theme
#' 
#' @description A custom ggplot2 theme for iHMP analysis visualizations.
#' Provides consistent styling across all plots in the package.
#' 
#' @return A ggplot2 theme object
#' @export
theme_ihmp <- function() {
  theme_minimal() +
    theme(
      # Text elements
      text = element_text(family = "Arial", size = 12),
      plot.title = element_text(size = 14, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 12, hjust = 0),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      
      # Panel elements
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.1),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.2),
      
      # Legend elements
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "grey50", linewidth = 0.2),
      legend.key = element_rect(fill = "white", color = NA),
      
      # Facet elements
      strip.background = element_rect(fill = "grey90", color = "grey50", linewidth = 0.2),
      strip.text = element_text(size = 10, face = "bold"),
      
      # Plot elements
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10)
    )
}

#' iHMP Color Palette
#' 
#' @description A custom color palette for iHMP analysis visualizations.
#' Provides consistent colors across all plots in the package.
#' 
#' @return A named vector of colors
#' @export
ihmp_colors <- function() {
  c(
    # Diagnosis colors
    "UC" = "#bb9df5",
    "nonIBD" = "#71d987",
    "CD" = "#ffc169",
    "IBD" = "#cb91ed",
    
    # Significance colors
    "significant" = "#e41a1c",
    "not_significant" = "grey70",
    
    # Mutation type colors
    "Missense" = "#1f77b4",
    "Silent" = "#2ca02c",
    "Intergenic" = "#ff7f0e",
    
    # Heatmap colors
    "heatmap_low" = "white",
    "heatmap_mid" = "#a8ddb5",
    "heatmap_high" = "#084081"
  )
} 
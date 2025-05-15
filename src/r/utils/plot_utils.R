#' Create Color Palette
#' 
#' @description Creates a color palette for plots.
#' 
#' @param n Number of colors to generate
#' @param palette Name of the color palette (default: "viridis")
#' 
#' @return Vector of colors
#' 
#' @export
create_palette <- function(n, palette = "viridis") {
  switch(palette,
         "viridis" = viridis::viridis(n),
         "magma" = viridis::magma(n),
         "plasma" = viridis::plasma(n),
         "inferno" = viridis::inferno(n),
         "cividis" = viridis::cividis(n),
         stop(sprintf("Unsupported palette: %s", palette)))
}

#' Format Plot
#' 
#' @description Formats a plot with consistent styling.
#' 
#' @param plot ggplot object
#' @param title Plot title
#' @param x_label X-axis label
#' @param y_label Y-axis label
#' @param theme Theme to use (default: theme_ihmp())
#' 
#' @return Formatted ggplot object
#' 
#' @export
format_plot <- function(plot, title = NULL, x_label = NULL, y_label = NULL,
                       theme = theme_ihmp()) {
  if (!inherits(plot, "ggplot")) {
    stop("plot must be a ggplot object")
  }
  
  plot <- plot + theme
  
  if (!is.null(title)) {
    plot <- plot + labs(title = title)
  }
  if (!is.null(x_label)) {
    plot <- plot + labs(x = x_label)
  }
  if (!is.null(y_label)) {
    plot <- plot + labs(y = y_label)
  }
  
  plot
}

#' Save Plot
#' 
#' @description Saves a plot to a file.
#' 
#' @param plot ggplot object
#' @param file_path Path to save the plot
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 8)
#' @param dpi Plot resolution (default: 300)
#' 
#' @export
save_plot <- function(plot, file_path, width = 12, height = 8, dpi = 300) {
  if (!inherits(plot, "ggplot")) {
    stop("plot must be a ggplot object")
  }
  
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  
  ggsave(file_path, plot, width = width, height = height, dpi = dpi)
}

#' Create Multi-Panel Plot
#' 
#' @description Creates a multi-panel plot from a list of plots.
#' 
#' @param plot_list List of ggplot objects
#' @param ncol Number of columns (default: 2)
#' @param title Overall plot title
#' 
#' @return Combined ggplot object
#' 
#' @export
create_multi_panel <- function(plot_list, ncol = 2, title = NULL) {
  if (!is.list(plot_list)) {
    stop("plot_list must be a list")
  }
  if (!all(sapply(plot_list, inherits, "ggplot"))) {
    stop("All elements in plot_list must be ggplot objects")
  }
  
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = ncol)
  
  if (!is.null(title)) {
    combined_plot <- combined_plot + 
      patchwork::plot_annotation(title = title)
  }
  
  combined_plot
} 
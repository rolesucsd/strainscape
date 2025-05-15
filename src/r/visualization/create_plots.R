#' Create Plots
#' 
#' @description Creates various plots for visualizing analysis results.
#' 
#' @param analysis_results List containing analysis results
#' @param output_dir Directory to save plots
#' @param plot_type Type of plot to create (default: "all")
#' @param theme Theme to use for plots (default: theme_ihmp())
#' 
#' @return A list of ggplot objects
#' 
#' @export
create_plots <- function(analysis_results, output_dir, plot_type = "all", theme = theme_ihmp()) {
  # Input validation
  if (!is.list(analysis_results)) {
    stop("analysis_results must be a list")
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  plots <- list()
  
  # Create coverage plots
  if (plot_type %in% c("all", "coverage")) {
    coverage_plot <- ggplot(analysis_results$coverage_stats,
                           aes(x = scaffold, y = mean_coverage)) +
      geom_bar(stat = "identity", fill = ihmp_colors()["coverage"]) +
      theme +
      labs(x = "Scaffold", y = "Mean Coverage",
           title = "Scaffold Coverage Distribution")
    
    plots$coverage <- coverage_plot
    ggsave(file.path(output_dir, "coverage_plot.pdf"),
           coverage_plot, width = 12, height = 8)
  }
  
  # Create mutation frequency plots
  if (plot_type %in% c("all", "mutations")) {
    mutation_plot <- ggplot(analysis_results$mutation_stats,
                           aes(x = mutation, y = mean_frequency)) +
      geom_point(aes(size = total_occurrences),
                 color = ihmp_colors()["mutation"]) +
      theme +
      labs(x = "Mutation", y = "Mean Frequency",
           title = "Mutation Frequency Distribution",
           size = "Total Occurrences")
    
    plots$mutation <- mutation_plot
    ggsave(file.path(output_dir, "mutation_plot.pdf"),
           mutation_plot, width = 12, height = 8)
  }
  
  # Create gene plots
  if (plot_type %in% c("all", "genes") && !is.null(analysis_results$gene_stats)) {
    gene_plot <- ggplot(analysis_results$gene_stats,
                        aes(x = gene, y = mean_frequency)) +
      geom_point(aes(size = total_mutations),
                 color = ihmp_colors()["gene"]) +
      theme +
      labs(x = "Gene", y = "Mean Frequency",
           title = "Gene Frequency Distribution",
           size = "Total Mutations")
    
    plots$gene <- gene_plot
    ggsave(file.path(output_dir, "gene_plot.pdf"),
           gene_plot, width = 12, height = 8)
  }
  
  # Create functional enrichment plots
  if (plot_type %in% c("all", "enrichment") && 
      !is.null(analysis_results$functional_enrichment)) {
    enrichment_plot <- ggplot(analysis_results$functional_enrichment,
                             aes(x = function, y = n_genes)) +
      geom_bar(stat = "identity", fill = ihmp_colors()["enrichment"]) +
      theme +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Function", y = "Number of Genes",
           title = "Functional Enrichment Analysis")
    
    plots$enrichment <- enrichment_plot
    ggsave(file.path(output_dir, "enrichment_plot.pdf"),
           enrichment_plot, width = 12, height = 8)
  }
  
  # Return plots
  plots
} 
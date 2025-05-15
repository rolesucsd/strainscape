#' Analyze Genes
#' 
#' @description Analyzes gene data from InStrain output to identify significant
#' changes in gene coverage and mutation patterns across timepoints.
#' 
#' @param gene_data A data frame containing gene information
#' @param metadata A data frame containing sample metadata
#' @param output_dir Directory to save output files
#' @param min_coverage Minimum coverage threshold (default: 5)
#' @param min_freq Minimum frequency threshold (default: 0.1)
#' @param p_threshold P-value threshold for significance (default: 0.05)
#' 
#' @return A list containing:
#' \itemize{
#'   \item significant_genes: Data frame of genes with significant changes
#'   \item gene_stats: Summary statistics of genes
#'   \item functional_enrichment: Functional enrichment analysis results
#' }
#' 
#' @export
analyze_genes <- function(gene_data, metadata, output_dir,
                         min_coverage = 5, min_freq = 0.1, p_threshold = 0.05) {
  # Input validation
  if (!is.data.frame(gene_data)) {
    stop("gene_data must be a data frame")
  }
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame")
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Filter genes by coverage and frequency
  filtered_genes <- gene_data %>%
    filter(coverage >= min_coverage, frequency >= min_freq)
  
  # Calculate gene statistics
  gene_stats <- filtered_genes %>%
    group_by(gene) %>%
    summarize(
      total_mutations = n(),
      mean_coverage = mean(coverage),
      mean_frequency = mean(frequency),
      max_frequency = max(frequency),
      min_frequency = min(frequency)
    )
  
  # Identify significant changes
  significant_genes <- filtered_genes %>%
    group_by(gene) %>%
    filter(n() > 1) %>%  # At least 2 timepoints
    summarize(
      coverage_change = max(coverage) - min(coverage),
      frequency_change = max(frequency) - min(frequency),
      p_value = t.test(frequency ~ timepoint)$p.value
    ) %>%
    filter(p_value < p_threshold)
  
  # Perform functional enrichment analysis
  if ("function" %in% colnames(gene_data)) {
    functional_enrichment <- significant_genes %>%
      left_join(gene_data %>% select(gene, function) %>% distinct(),
                by = "gene") %>%
      group_by(function) %>%
      summarize(
        n_genes = n(),
        mean_frequency_change = mean(frequency_change),
        mean_coverage_change = mean(coverage_change)
      ) %>%
      arrange(desc(n_genes))
    
    # Save functional enrichment results
    write.csv(functional_enrichment,
              file.path(output_dir, "functional_enrichment.csv"),
              row.names = FALSE)
  } else {
    functional_enrichment <- NULL
  }
  
  # Save results
  write.csv(significant_genes,
            file.path(output_dir, "significant_genes.csv"),
            row.names = FALSE)
  write.csv(gene_stats,
            file.path(output_dir, "gene_stats.csv"),
            row.names = FALSE)
  
  # Return results
  list(
    significant_genes = significant_genes,
    gene_stats = gene_stats,
    functional_enrichment = functional_enrichment
  )
} 
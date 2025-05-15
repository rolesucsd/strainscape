#' Analyze Scaffolds
#' 
#' @description Analyzes scaffold data from InStrain output to identify significant changes
#' in scaffold coverage and mutation patterns across timepoints.
#' 
#' @param scaffold_data A data frame containing scaffold information
#' @param metadata A data frame containing sample metadata
#' @param output_dir Directory to save output files
#' @param min_coverage Minimum coverage threshold (default: 5)
#' @param min_freq Minimum frequency threshold (default: 0.1)
#' 
#' @return A list containing:
#' \itemize{
#'   \item significant_scaffolds: Data frame of scaffolds with significant changes
#'   \item coverage_stats: Summary statistics of scaffold coverage
#'   \item mutation_stats: Summary statistics of mutations
#' }
#' 
#' @export
analyze_scaffolds <- function(scaffold_data, metadata, output_dir, 
                            min_coverage = 5, min_freq = 0.1) {
  # Input validation
  if (!is.data.frame(scaffold_data)) {
    stop("scaffold_data must be a data frame")
  }
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame")
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Filter scaffolds by coverage and frequency
  filtered_scaffolds <- scaffold_data %>%
    filter(coverage >= min_coverage, frequency >= min_freq)
  
  # Calculate statistics
  coverage_stats <- filtered_scaffolds %>%
    group_by(scaffold) %>%
    summarize(
      mean_coverage = mean(coverage),
      sd_coverage = sd(coverage),
      max_coverage = max(coverage),
      min_coverage = min(coverage)
    )
  
  mutation_stats <- filtered_scaffolds %>%
    group_by(scaffold) %>%
    summarize(
      total_mutations = n(),
      unique_mutations = n_distinct(mutation),
      mean_frequency = mean(frequency)
    )
  
  # Identify significant changes
  significant_scaffolds <- filtered_scaffolds %>%
    group_by(scaffold) %>%
    filter(n() > 1) %>%  # At least 2 timepoints
    summarize(
      coverage_change = max(coverage) - min(coverage),
      mutation_change = max(frequency) - min(frequency)
    ) %>%
    filter(coverage_change > 0 | mutation_change > 0)
  
  # Save results
  write.csv(significant_scaffolds, 
            file.path(output_dir, "significant_scaffolds.csv"),
            row.names = FALSE)
  write.csv(coverage_stats,
            file.path(output_dir, "coverage_stats.csv"),
            row.names = FALSE)
  write.csv(mutation_stats,
            file.path(output_dir, "mutation_stats.csv"),
            row.names = FALSE)
  
  # Return results
  list(
    significant_scaffolds = significant_scaffolds,
    coverage_stats = coverage_stats,
    mutation_stats = mutation_stats
  )
} 
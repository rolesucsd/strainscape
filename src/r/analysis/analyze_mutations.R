#' Analyze Mutations
#' 
#' @description Analyzes mutation data from InStrain output to identify significant
#' changes in mutation patterns across timepoints.
#' 
#' @param mutation_data A data frame containing mutation information
#' @param metadata A data frame containing sample metadata
#' @param output_dir Directory to save output files
#' @param min_coverage Minimum coverage threshold (default: 5)
#' @param min_freq Minimum frequency threshold (default: 0.1)
#' @param p_threshold P-value threshold for significance (default: 0.05)
#' 
#' @return A list containing:
#' \itemize{
#'   \item significant_mutations: Data frame of mutations with significant changes
#'   \item mutation_stats: Summary statistics of mutations
#'   \item timepoint_comparisons: Pairwise comparisons between timepoints
#' }
#' 
#' @export
analyze_mutations <- function(mutation_data, metadata, output_dir,
                            min_coverage = 5, min_freq = 0.1, p_threshold = 0.05) {
  # Input validation
  if (!is.data.frame(mutation_data)) {
    stop("mutation_data must be a data frame")
  }
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame")
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Filter mutations by coverage and frequency
  filtered_mutations <- mutation_data %>%
    filter(coverage >= min_coverage, frequency >= min_freq)
  
  # Calculate mutation statistics
  mutation_stats <- filtered_mutations %>%
    group_by(mutation) %>%
    summarize(
      total_occurrences = n(),
      mean_frequency = mean(frequency),
      max_frequency = max(frequency),
      min_frequency = min(frequency)
    )
  
  # Identify significant changes
  significant_mutations <- filtered_mutations %>%
    group_by(mutation) %>%
    filter(n() > 1) %>%  # At least 2 timepoints
    summarize(
      frequency_change = max(frequency) - min(frequency),
      p_value = t.test(frequency ~ timepoint)$p.value
    ) %>%
    filter(p_value < p_threshold)
  
  # Perform pairwise comparisons
  timepoint_comparisons <- filtered_mutations %>%
    group_by(mutation) %>%
    filter(n() > 1) %>%
    do({
      timepoints <- unique(.$timepoint)
      comparisons <- combn(timepoints, 2, simplify = FALSE)
      results <- lapply(comparisons, function(tp) {
        data.frame(
          mutation = .$mutation[1],
          timepoint1 = tp[1],
          timepoint2 = tp[2],
          freq_diff = diff(.$frequency[.$timepoint %in% tp]),
          p_value = t.test(.$frequency[.$timepoint == tp[1]],
                          .$frequency[.$timepoint == tp[2]])$p.value
        )
      })
      do.call(rbind, results)
    })
  
  # Save results
  write.csv(significant_mutations,
            file.path(output_dir, "significant_mutations.csv"),
            row.names = FALSE)
  write.csv(mutation_stats,
            file.path(output_dir, "mutation_stats.csv"),
            row.names = FALSE)
  write.csv(timepoint_comparisons,
            file.path(output_dir, "timepoint_comparisons.csv"),
            row.names = FALSE)
  
  # Return results
  list(
    significant_mutations = significant_mutations,
    mutation_stats = mutation_stats,
    timepoint_comparisons = timepoint_comparisons
  )
} 
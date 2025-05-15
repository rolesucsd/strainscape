# -----------------------------
# iHMP Integrated Microbial Evolution Analysis
# -----------------------------
# This script integrates multiple analyses to understand microbial evolution in IBD:
# 1. Combines dN/dS, nucleotide diversity, and mutation accrual patterns
# 2. Analyzes selection pressures and adaptation patterns
# 3. Correlates evolutionary patterns with clinical outcomes
# 4. Identifies potential biomarkers of disease progression

# Libraries
library(tidyverse)
library(data.table)    # for fast I/O and operations
library(progress)      # for progress bars
library(ggpubr)        # for statistical tests
library(ggdist)        # for raincloud plots
library(phyloseq)      # for microbial community analysis
library(vegan)         # for diversity metrics
library(corrr)         # for correlation analysis

# ----- Configuration & Paths -----
ROOT <- "../iHMP"
FIG_DIR <- file.path(ROOT, "samples", "wolr2", "figures", "integrated_analysis")

# ----- Constants with Biological Significance -----
# MIN_COVERAGE: Minimum coverage for reliable variant calling (20x)
# - Ensures accurate detection of mutations
# - Reduces false positive variant calls
MIN_COVERAGE <- 20

# MIN_ABUNDANCE: Minimum relative abundance for species analysis (0.1%)
# - Ensures sufficient representation for evolutionary analysis
# - Reduces noise from low-abundance species
MIN_ABUNDANCE <- 0.001

# ----- Helper Functions -----
integrate_evolutionary_metrics <- function(dnds_data, diversity_data, mutation_data) {
  # Integrate different evolutionary metrics with progress bar
  pb <- progress_bar$new(total = 3, format = "Integrating metrics [:bar] :percent")
  pb$tick()
  
  # Combine dN/dS and diversity data
  integrated_data <- dnds_data[
    diversity_data$species_df,
    on = c("Species", "Sample" = "Participant ID")
  ]
  
  pb$tick()
  
  # Add mutation accrual patterns
  integrated_data <- integrated_data[
    mutation_data,
    on = c("Species", "Sample")
  ]
  
  pb$tick()
  
  integrated_data
}

analyze_selection_pressures <- function(integrated_data) {
  # Analyze selection pressures with progress bar
  pb <- progress_bar$new(total = 3, format = "Analyzing selection [:bar] :percent")
  pb$tick()
  
  # Calculate selection metrics
  selection_metrics <- integrated_data[
    , .(
      # Selection pressure metrics
      dnds_ratio = mean(fraction_dn / fraction_ds, na.rm = TRUE),
      diversity_trend = coef(lm(mean_div ~ week_num))[2],
      mutation_rate = mean(abs_rate_per_wk, na.rm = TRUE),
      
      # Adaptation metrics
      adaptive_sweeps = sum(OLS_pvalue <= 0.05 & abs(OLS_slope) > 0.001, na.rm = TRUE),
      selective_sweeps = sum(delta > 0.8, na.rm = TRUE),
      
      # Population metrics
      effective_pop_size = 1 / (4 * mean(nucl_diversity, na.rm = TRUE)),
      population_bottleneck = min(nucl_diversity, na.rm = TRUE) / max(nucl_diversity, na.rm = TRUE)
    ),
    by = .(Species, Sample, diagnosis)
  ]
  
  pb$tick()
  
  # Add statistical significance
  selection_metrics <- selection_metrics[
    , `:=`(
      selection_significant = dnds_ratio > 1 & adaptive_sweeps > 0,
      adaptation_significant = selective_sweeps > 0 & diversity_trend < 0
    )
  ]
  
  pb$tick()
  
  selection_metrics
}

correlate_with_clinical_outcomes <- function(selection_metrics, clinical_data) {
  # Correlate evolutionary patterns with clinical outcomes
  pb <- progress_bar$new(total = 2, format = "Correlating with outcomes [:bar] :percent")
  pb$tick()
  
  # Join with clinical data
  clinical_correlations <- selection_metrics[
    clinical_data,
    on = "Sample"
  ]
  
  # Calculate correlations
  correlations <- clinical_correlations[
    , .(
      correlation = cor(dnds_ratio, disease_score, method = "spearman"),
      p_value = cor.test(dnds_ratio, disease_score, method = "spearman")$p.value
    ),
    by = .(Species, diagnosis)
  ]
  
  pb$tick()
  
  correlations
}

create_integrated_visualizations <- function(selection_metrics, correlations, fig_dir) {
  # Create integrated visualizations
  pb <- progress_bar$new(total = 3, format = "Creating visualizations [:bar] :percent")
  
  # 1. Selection pressure heatmap
  pb$tick()
  p1 <- selection_metrics[
    , .(mean_dnds = mean(dnds_ratio, na.rm = TRUE)),
    by = .(Species, diagnosis)
  ] %>%
    ggplot(aes(diagnosis, Species, fill = mean_dnds)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 1,
      name = "Mean dN/dS"
    ) +
    labs(title = "Selection Pressures Across Species",
         subtitle = "dN/dS ratios by diagnosis and species",
         x = NULL,
         y = NULL) +
    custom_theme
  
  # 2. Adaptation patterns
  pb$tick()
  p2 <- selection_metrics[
    , .(
      n_adaptive = sum(selection_significant),
      n_adapted = sum(adaptation_significant)
    ),
    by = .(Species, diagnosis)
  ] %>%
    melt(
      measure.vars = c("n_adaptive", "n_adapted"),
      variable.name = "metric",
      value.name = "count"
    ) %>%
    ggplot(aes(diagnosis, count, fill = metric)) +
    geom_col(position = "dodge") +
    facet_wrap(~ Species) +
    labs(title = "Adaptation Patterns",
         subtitle = "Number of adaptive and adapted strains by species",
         x = NULL,
         y = "Count") +
    custom_theme
  
  # 3. Clinical correlations
  pb$tick()
  p3 <- correlations[
    p_value < 0.05
  ] %>%
    ggplot(aes(Species, correlation, fill = diagnosis)) +
    geom_col(position = "dodge") +
    coord_flip() +
    labs(title = "Clinical Correlations",
         subtitle = "Significant correlations with disease severity",
         x = NULL,
         y = "Spearman correlation") +
    custom_theme
  
  # Save plots
  ggsave(file.path(fig_dir, "selection_pressures.png"), p1, width = 10, height = 8)
  ggsave(file.path(fig_dir, "adaptation_patterns.png"), p2, width = 12, height = 10)
  ggsave(file.path(fig_dir, "clinical_correlations.png"), p3, width = 8, height = 6)
}

# ----- Main Analysis -----
# Integrate data from previous analyses
integrated_data <- integrate_evolutionary_metrics(dnds_ratios, diversity_data, snp_metrics)

# Analyze selection pressures
selection_metrics <- analyze_selection_pressures(integrated_data)

# Correlate with clinical outcomes
clinical_correlations <- correlate_with_clinical_outcomes(selection_metrics, iHMP_metadata_disease)

# Create visualizations
create_integrated_visualizations(selection_metrics, clinical_correlations, FIG_DIR)

# ----- Output Summary Statistics -----
# Write summary statistics
fwrite(selection_metrics, file.path(FIG_DIR, "selection_metrics.txt"), sep = "\t")
fwrite(clinical_correlations, file.path(FIG_DIR, "clinical_correlations.txt"), sep = "\t")

# ----- Additional Analysis Suggestions -----
# 1. Time series analysis of selection pressures
# 2. Phylogenetic analysis of adapted strains
# 3. Functional enrichment of selected mutations
# 4. Network analysis of co-evolving species
# 5. Machine learning for prediction of disease progression 
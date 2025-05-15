# -----------------------------
# iHMP dN/dS Analysis
# -----------------------------
# This script analyzes the ratio of non-synonymous to synonymous mutations (dN/dS):
# 1. Calculates dN/dS ratios per strain and species
# 2. Analyzes mutation patterns across different diagnoses
# 3. Generates visualizations of dN/dS distributions
# 4. Identifies potential selection pressures in different conditions

# Libraries
library(tidyverse)
library(data.table)    # for fast I/O and operations
library(progress)      # for progress bars
library(ggdist)        # for raincloud plots

# ----- Configuration & Paths -----
FIG_DIR <- "../iHMP/samples/wolr2/figures/dnds"

# ----- Constants with Biological Significance -----
# MIN_MUTATIONS: Minimum number of mutations required for reliable dN/dS calculation
# - Ensures statistical significance of the ratio
# - Reduces noise from small sample sizes
MIN_MUTATIONS <- 10

# ----- Helper Functions -----
calculate_dnds_ratios <- function(combined_df_subset_metadata_summary) {
  # Calculate dN/dS ratios with progress bar
  pb <- progress_bar$new(total = 3, format = "Calculating dN/dS ratios [:bar] :percent")
  pb$tick()
  
  # Filter and prepare data
  ds_dn <- combined_df_subset_metadata_summary[
    coding == "genic",
    .(Species, Sample, Chromosome, Position, Mutation_Type)
  ] %>%
    unique() %>%
    .[, category := fcase(
      Mutation_Type == "Silent", "ds",
      Mutation_Type %in% c("Missense", "Nonsense"), "dn",
      default = NA_character_
    )] %>%
    .[!is.na(category)]
  
  pb$tick()
  
  # Calculate ratios
  ds_dn <- ds_dn[
    , .N,
    by = .(Species, Sample, category)
  ] %>%
    dcast(
      Species + Sample ~ category,
      value.var = "N",
      fill = 0
    ) %>%
    .[, `:=`(
      total = dn + ds,
      fraction_dn = dn / total,
      fraction_ds = ds / total
    )]
  
  pb$tick()
  
  # Join with metadata
  ds_dn <- ds_dn[
    iHMP_metadata_disease_summarized,
    on = "Sample"
  ]
  
  ds_dn
}

create_dnds_visualizations <- function(ds_dn, fig_dir) {
  # Create visualizations with progress bar
  pb <- progress_bar$new(total = 2, format = "Creating visualizations [:bar] :percent")
  
  # Define colors
  colors <- c("UC"="#bb9df5", "nonIBD"="#71d987", "CD"="#ffc169", "IBD"="#cb91ed")
  
  # 1. Species-specific boxplots
  pb$tick()
  p1 <- ggplot(ds_dn[total >= MIN_MUTATIONS],
               aes(x = diagnosis, y = fraction_dn, fill = diagnosis)) +
    geom_boxplot() +
    geom_jitter(width = 0.05) +
    facet_wrap(~Species) +
    scale_fill_manual(values = colors) +
    labs(title = "dN/dS Ratio Distribution by Species",
         subtitle = "Showing only strains with ≥10 mutations",
         x = NULL,
         y = "dN/dS ratio") +
    custom_theme
  
  ggsave(file.path(fig_dir, "dnds_species.png"),
         p1, height = 15, width = 15, units = "in")
  
  # 2. Overall distribution
  pb$tick()
  p2 <- ggplot(ds_dn[diagnosis %in% c("nonIBD", "CD", "UC")],
               aes(x = diagnosis, y = fraction_dn, fill = diagnosis, color = diagnosis)) +
    geom_jitter(width = .05, size = 1, alpha = .75) +
    ggdist::stat_halfeye(adjust = .5, width = .5, .width = 0, alpha = .8, color = "black") +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    labs(title = "Overall dN/dS Ratio Distribution",
         subtitle = "By diagnosis group",
         x = NULL,
         y = "dN/dS ratio") +
    custom_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(fig_dir, "dnds.png"),
         p2, height = 3.5, width = 2, units = "in")
}

# ----- Main Analysis -----
# Calculate dN/dS ratios
ds_dn_ratios <- calculate_dnds_ratios(combined_df_subset_metadata_summary)

# Create visualizations
create_dnds_visualizations(ds_dn_ratios, FIG_DIR)

# ----- Optional: Additional Analysis -----
# Add any additional analysis here, such as:
# - Statistical tests for differences between groups
# - Correlation with other metrics
# - Time series analysis of dN/dS ratios

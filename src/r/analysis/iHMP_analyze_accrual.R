#' Analyze mutation accrual patterns in iHMP dataset
#' 
#' @description This function analyzes mutation accrual patterns in the iHMP dataset,
#' calculating mutation metrics, adaptive sweep rates, and generating visualizations.
#' 
#' @param combined_df Data frame containing combined mutation data
#' @param metadata Data frame containing sample metadata
#' @param outdir Directory to save output files
#' @param low_thr Threshold for low frequency mutations (default: 0.10)
#' @param high_thr Threshold for high frequency mutations (default: 0.80)
#' @param sweep_delta Threshold for adaptive sweep detection (default: 0.30)
#' 
#' @return List containing processed data and summary metrics
#' @export
analyze_mutation_accrual <- function(combined_df, metadata, outdir,
                                   low_thr = 0.10, high_thr = 0.80, sweep_delta = 0.30) {
  # Input validation
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  # Process mutation data
  snp_metrics <- process_mutation_data(combined_df)
  
  # Calculate summary metrics
  summary_metrics <- calculate_summary_metrics(snp_metrics)
  
  # Create visualizations
  create_visualizations(snp_metrics, summary_metrics, outdir)
  
  # Return results
  list(
    snp_metrics = snp_metrics,
    summary_metrics = summary_metrics
  )
}

#' Process mutation data
#' @keywords internal
process_mutation_data <- function(combined_df) {
  # Process mutation data with progress bar
  pb <- progress_bar$new(total = 3, format = "Processing mutation data [:bar] :percent")
  pb$tick()
  
  # Clean and prepare data
  df <- combined_df[!is.na(percentage)] %>%
    .[, `:=`(
      week = as.numeric(week),
      percentage = as.numeric(percentage)
    )]
  
  pb$tick()
  
  # Calculate SNP metrics
  snp_metrics <- df[
    order(Sample, Species, snp_id, week),
    .(
      max_prop = max(percentage, na.rm = TRUE),
      week_max_prop = week[which.max(percentage)],
      min_prop = min(percentage, na.rm = TRUE),
      week_min_prop = week[which.min(percentage)],
      delta = max(percentage, na.rm = TRUE) - min(percentage, na.rm = TRUE)
    ),
    by = .(Sample, Species, snp_id)
  ]
  
  pb$tick()
  
  # Join metadata and clean mutation types
  snp_metrics <- snp_metrics[combined_df_subset_metadata_summary, on = c("Sample", "Species")] %>%
    .[!is.na(OLS_slope)] %>%
    .[, `:=`(
      Mutation_Type = fcase(
        Mutation_Type == "intergenic", "Intergenic",
        Mutation_Type == "Nonsense", "Missense",
        Mutation_Type %like% "Position", "Silent",
        default = Mutation_Type
      ),
      OLS_slope = as.numeric(OLS_slope),
      Species = sub(".*s__", "", Species),
      diagnosis = factor(diagnosis, levels = c("nonIBD", "CD", "UC"))
    )]
  
  snp_metrics
}

#' Calculate summary metrics
#' @keywords internal
calculate_summary_metrics <- function(snp_metrics) {
  # Calculate summary metrics with progress bar
  pb <- progress_bar$new(total = 2, format = "Calculating summary metrics [:bar] :percent")
  pb$tick()
  
  # Calculate summary statistics
  summary_metrics <- snp_metrics[
    , .(
      new_snv = sum(delta > 0.8, na.rm = TRUE),
      time_window_wk = diff(range(c(week_min_prop, week_max_prop))),
      abs_rate_per_wk = sum(delta > 0.8, na.rm = TRUE) / diff(range(c(week_min_prop, week_max_prop))),
      n_adaptive = sum(OLS_pvalue <= 0.10, na.rm = TRUE),
      mean_change_freq = mean(delta, na.rm = TRUE),
      mean_slope = mean(as.numeric(OLS_slope), na.rm = TRUE)
    ),
    by = .(Sample, Species, diagnosis)
  ]
  
  pb$tick()
  
  # Filter out extreme values
  summary_metrics <- summary_metrics[new_snv < 500]
  
  summary_metrics
}

#' Create visualizations
#' @keywords internal
create_visualizations <- function(snp_metrics, snp_metrics_summary, outdir) {
  # Create visualizations with progress bar
  pb <- progress_bar$new(total = 5, format = "Creating visualizations [:bar] :percent")
  
  # Define colors
  colors <- c("UC"="#bb9df5", "nonIBD"="#71d987", "CD"="#ffc169", "IBD"="#cb91ed")
  
  # 1. Raincloud plot of sweep rates
  pb$tick()
  p1 <- ggplot(snp_metrics_summary, 
               aes(diagnosis, mean_change_freq, colour = diagnosis, fill = diagnosis)) +
    ggdist::stat_halfeye(adjust = .5, width = .5, .width = 0, alpha = .8, color = "black") +
    geom_jitter(width = .05, size = 1, alpha = .75) +
    scale_fill_manual(values = colors) +
    scale_colour_manual(values = colors) +
    labs(title = "Distribution of Adaptive Sweep Rates",
         subtitle = "By diagnosis group",
         y = expression(log[10] * " adaptive-sweep rate (wk"^-1*")"),
         x = NULL) +
    theme_ihmp() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(outdir, "Fig_Raincloud_SweepRate.png"),
         p1, width = 2, height = 3.5, units = "in", dpi=600)
  
  # 2. Heatmap of species rates
  pb$tick()
  species_counts <- snp_metrics_summary[, .N, by = Species][N >= 3]
  species_pat_counts <- snp_metrics_summary[
    Species %in% species_counts$Species,
    .(n_patients = uniqueN(Sample)),
    by = Species
  ]
  
  heat_df <- snp_metrics_summary[
    Species %in% species_pat_counts$Species,
    .(mean_rate = mean(abs_rate_per_wk, na.rm = TRUE)),
    by = .(Species, diagnosis)
  ]
  
  mat <- dcast(heat_df, Species ~ diagnosis, value.var = "mean_rate") %>%
    as.matrix()
  rownames(mat) <- mat[,1]
  mat <- mat[,-1]
  mat[is.na(mat)] <- 0
  mat <- log10(mat + 1e-3)
  
  ha <- rowAnnotation(
    Patients = species_pat_counts$n_patients[match(rownames(mat), species_pat_counts$Species)],
    col = list(Patients = colorRamp2(c(1, max(species_pat_counts$n_patients)),
                                   c("white","black")))
  )
  
  col_fun <- colorRamp2(
    c(-3, -2, -1, 0, 1),
    c("white", "#e0f3db", "#a8ddb5", "#43a2ca", "#084081")
  )
  
  png(file.path(outdir, "Fig_Heat_SpeciesRate.png"),
      width = 5, height = 6, units = "in", res=600)
  Heatmap(mat,
          name = "-log10 sweep rate",
          col = col_fun,
          right_annotation = ha,
          cluster_rows = TRUE,
          cluster_columns = FALSE)
  dev.off()
  
  # 3. Mutation type proportions
  pb$tick()
  mut_prop <- snp_metrics[
    Mutation_Type %in% c("Missense","Intergenic","Silent"),
    .N,
    by = .(diagnosis, Mutation_Type)
  ][, prop := N/sum(N), by = diagnosis]
  
  p3 <- ggplot(mut_prop, aes(diagnosis, prop, fill = Mutation_Type)) +
    geom_col(width = .6) +
    scale_fill_brewer(palette="Set2") +
    labs(title = "Distribution of Mutation Types",
         subtitle = "By diagnosis group",
         y = "Proportion of SNPs",
         x = NULL,
         fill = NULL) +
    theme_ihmp() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(outdir, "Fig_Metrics.png"),
         p3, width = 2, height = 3.5, units = "in", dpi=600)
  
  # 4. Disease score correlation
  pb$tick()
  inflam_long <- snp_metrics_summary[
    iHMP_metadata_disease[, .(Sample = `Participant ID`, week = week_num, hbi, sccai)],
    on = "Sample"
  ] %>%
    melt(
      measure.vars = c("hbi", "sccai"),
      variable.name = "score_type",
      value.name = "disease_score"
    ) %>%
    .[, .(
      disease_score = mean(disease_score, na.rm = TRUE),
      abs_rate_per_wk = mean(abs_rate_per_wk, na.rm = TRUE)
    ),
      by = .(Sample, diagnosis, score_type)]
  
  p4 <- ggplot(inflam_long, 
               aes(x = disease_score, y = abs_rate_per_wk, colour = diagnosis)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
    scale_colour_manual(values = colors) +
    scale_y_log10() +
    labs(title = "Correlation with Disease Scores",
         subtitle = "Adaptive sweep rates vs. HBI/SCCAI scores",
         x = "Disease Score (HBI or SCCAI)",
         y = "Adaptive-sweep rate/wk") +
    theme_ihmp() +
    ggpubr::stat_cor(method = "spearman", label.y.npc = 0.95)
  
  # 5. Slope vs p-value plot
  pb$tick()
  p5 <- snp_metrics[
    , `:=`(
      log_slope = log2(abs(OLS_slope*100)+1e-4),
      neglog_p = -log10(OLS_pvalue)
    )
  ] %>%
    ggplot(aes(log_slope, neglog_p, colour = diagnosis)) +
    geom_point(size = .4, alpha = .6) +
    scale_colour_manual(values = colors) +
    facet_wrap(~ Species) +
    labs(title = "Mutation Significance vs. Rate of Change",
         subtitle = "By species and diagnosis",
         x = expression(log[2] * "|" * slope * "|"),
         y = expression(-log[10]~italic(p))) +
    theme_ihmp()
} 
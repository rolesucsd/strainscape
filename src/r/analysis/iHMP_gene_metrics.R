#' Calculate gene-level metrics from SNP data
#' 
#' @description This function calculates gene-level metrics from SNP data,
#' including mutation rates, slopes, and statistical significance.
#' 
#' @param snp_metrics Data frame containing SNP-level metrics
#' @param outdir Directory to save output files
#' @param p_thresh P-value threshold for significance (default: 0.1)
#' @param enrich_thresh Enrichment threshold (default: 0.2)
#' 
#' @return Data frame containing gene-level metrics and statistics
#' @export
calculate_gene_metrics <- function(snp_metrics, outdir, p_thresh = 0.1, enrich_thresh = 0.2) {
  # Input validation
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  # Create gene_metrics dataframe
  gene_metrics <- snp_metrics %>%
    mutate(disease_group = if_else(diagnosis %in% c("UC", "CD"), "IBD", "nonIBD")) %>%
    group_by(Sample, Species, disease_group, KEGG, Mutation_Type) %>%
    summarise(
      n               = n(),                              # true SNP count
      mean_delta      = mean(delta,          na.rm = TRUE),
      mean_slope      = mean(abs(OLS_slope), na.rm = TRUE),
      .groups         = "drop"
    )
  
  # Summarize data by disease category
  gene_disease <- gene_metrics %>%
    group_by(KEGG, disease_group, Mutation_Type) %>%
    summarise(
      mean_delta = mean(mean_delta, na.rm = TRUE),
      mean_slope = mean(mean_slope, na.rm = TRUE),
      n          = sum(n),           # total SNPs in that gene/category
      .groups    = "drop"
    )
  
  # Calculate ranked metrics
  ranked_all <- calculate_ranked_metrics(gene_disease, gene_metrics)
  
  # Create volcano plot
  create_volcano_plot(ranked_all, outdir, p_thresh, enrich_thresh)
  
  # Return results
  ranked_all
}

#' Calculate ranked metrics for genes
#' @keywords internal
calculate_ranked_metrics <- function(gene_disease, gene_metrics) {
  # Helper: robust z-score
  safe_z <- function(x) {
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    (x - m) / s
  }
  
  # Core routine: returns a ranked table for ONE mutation class
  rank_one_class <- function(df_class) {
    df_class %>%
      count(KEGG, name = "n_groups") %>%
      filter(n_groups == 2) %>%
      inner_join(df_class, by = "KEGG") %>%               # only genes present in both groups
      pivot_wider(
        names_from  = disease_group,
        values_from = c(mean_delta, mean_slope, n),
        names_glue  = "{.value}_{disease_group}"
      ) %>%
      mutate(
        z_delta = safe_z(mean_delta_IBD - mean_delta_nonIBD),
        z_slope = safe_z(mean_slope_IBD - mean_slope_nonIBD),
        z_total = z_delta + z_slope
      ) %>%
      arrange(desc(z_total))
  }
  
  # Calculate ranked metrics for all mutation types
  ranked_all <- gene_disease %>%                          # already summarised
    group_by(Mutation_Type) %>%
    group_modify(~ rank_one_class(.x)) %>%                # .x is the subgroup
    ungroup()
  
  # Filter for genes with sufficient data
  keepers <- gene_metrics %>%
    group_by(KEGG, Mutation_Type, disease_group) %>%
    summarise(n = sum(!is.na(mean_delta)), .groups = "drop") %>%
    pivot_wider(names_from = disease_group, values_from = n, values_fill = 0) %>%
    filter(IBD >= 2, nonIBD >= 2) %>%
    select(KEGG, Mutation_Type)
  
  # Calculate statistical significance
  safe_wilcox <- function(x, y) {
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    if(length(x) >= 2 && length(y) >= 2) {
      wilcox.test(x, y, exact = FALSE)$p.value
    } else {
      NA_real_                                             # mark insufficient data
    }
  }
  
  # Add statistical significance
  ranked_all <- ranked_all %>%
    semi_join(keepers, by = c("KEGG", "Mutation_Type")) %>%
    rowwise() %>%
    mutate(
      p_val = safe_wilcox(
        gene_metrics$mean_delta[
          gene_metrics$KEGG           == KEGG &
            gene_metrics$disease_group  == "IBD" &
            gene_metrics$Mutation_Type  == Mutation_Type
        ],
        gene_metrics$mean_delta[
          gene_metrics$KEGG           == KEGG &
            gene_metrics$disease_group  == "nonIBD" &
            gene_metrics$Mutation_Type  == Mutation_Type
        ]
      )
    ) %>%
    ungroup() %>%
    mutate(
      enrich = log10((mean_delta_IBD + 1e-4) / (mean_delta_nonIBD + 1e-4)),
      log10P = -log10(p_val)
    )
  
  ranked_all
}

#' Create volcano plot of gene metrics
#' @keywords internal
create_volcano_plot <- function(ranked_all, outdir, p_thresh, enrich_thresh) {
  # Add significance flag
  ranked_all <- ranked_all %>%
    mutate(sig = p_val < p_thresh & abs(enrich) > enrich_thresh)
  
  # Define colors
  colors_scatter <- c("FALSE" = "grey70", "TRUE" = "red")
  
  # Create plot
  p <- ggplot(ranked_all[ranked_all$Mutation_Type == "Missense",],
              aes(x = enrich, y = log10P, colour = sig)) +
    geom_jitter(size = 3, height= 0.02) +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed") +
    geom_vline(xintercept = c(-enrich_thresh, enrich_thresh), linetype = "dashed") +
    scale_colour_manual(values = colors_scatter, guide = FALSE) +
    geom_text_repel(
      data = subset(ranked_all[ranked_all$Mutation_Type == "Missense",], sig),
      aes(label = KEGG), size = 3, max.overlaps = Inf
    ) +
    facet_wrap(~ Mutation_Type, scales = "free_x") +
    labs(x = "log10(enrichment IBD / nonIBD)",
         y = "-log10(p value)") +
    theme_ihmp()
  
  # Save plot
  ggsave(file.path(outdir, "Gene_Metrics_Volcano.png"),
         p, width = 4, height = 3.5, units = "in", dpi=600)
} 
# -----------------------------
# Iterative Sample × Species Mutation Analysis (Top 20 by Significance)
# -----------------------------
# Libraries
library(tidyverse)
library(stringr)
library(ggrepel)
library(RColorBrewer)  # for qualitative palettes


# ----- Configuration & Paths -----
# Data frames assumed already loaded:
#   nucleotide_div_filter: df with Sample (`Participant ID`) and Species
#   combined_df_subset: mutation metrics
#   scaffold_data_edit: metadata (if needed)
fig_dir <- "../iHMP/samples/wolr2/figures/isolates"

# ----- Thresholds -----
MIN_SLOPE <- 0.01   # absolute slope cutoff
MAX_PVAL  <- 0.05   # p-value significance cutoff
MAX_LINES <- 20     # max lines per plot

# ----- Prepare combinations -----
combos <- nucleotide_div_filter %>%
  distinct(`Participant ID`, Species) %>%
  arrange(`Participant ID`, Species)

# ----- Iterate over each Sample–Species combo -----
for(i in seq_len(nrow(combos))) {
  samp <- combos$`Participant ID`[i]
  sp   <- combos$Species[i]
  message("Processing Sample=", samp, " Species=", sp)
  
  # Subset relevant mutations
  df <- combined_df_subset %>%
    filter(Sample == samp,
           str_detect(Species, fixed(sp)))
  if(nrow(df) == 0) next
  
  # Annotate Mutation_Type and unified Product
  df <- df %>%
    mutate(
      Mutation_Type = if_else(is.na(Mutation_Type) | Mutation_Type == "", 
                              "Intergenic", Mutation_Type),
      Product = case_when(
        Mutation_Type %in% c("Missense", "Silent") ~ Matched_Product,
        Mutation_Type == "Intergenic" ~ paste(upstream_Product, downstream_Product, sep = "/"),
        TRUE ~ Matched_Product
      )
    )
  
  # Convert types and adjust percentages
  df <- df %>%
    mutate(
      week       = as.numeric(week),
      percentage = as.numeric(percentage),
      OLS_slope  = as.numeric(OLS_slope),
      OLS_pvalue = as.numeric(OLS_pvalue)
    ) %>%
    group_by(Position) %>%
    mutate(
      percentage = if_else(OLS_slope < 0, 1 - percentage, percentage)
    ) %>%
    ungroup()
  
  # Filter by thresholds
  df <- df %>%
    filter(abs(OLS_slope) > MIN_SLOPE,
           OLS_pvalue <= MAX_PVAL)
  if(nrow(df) == 0) next
  
  # Choose best Position per Product (highest |slope|), keep all weeks
  best_positions <- df %>%
    distinct(Product, Position, OLS_slope) %>%
    mutate(abs_slope = abs(OLS_slope)) %>%
    group_by(Product) %>%
    slice_max(abs_slope, n = 1, with_ties = FALSE) %>%
    pull(Position)
  df <- df %>% filter(Position %in% best_positions)
  if(nrow(df) == 0) next
  
  # Limit to top MIN(MAX_LINES, distinct Positions) by smallest p-value
  top_positions <- df %>%
    distinct(Position, OLS_pvalue) %>%
    arrange(OLS_pvalue) %>%
    slice_head(n = MAX_LINES) %>%
    pull(Position)
  cutoff_p <- max(df$OLS_pvalue[df$Position %in% top_positions])
  df <- df %>% filter(Position %in% top_positions)
  if(nrow(df) == 0) next
  
  # Determine label positions at median week per Position, keeping facet context
  label_df <- df %>%
    group_by(Mutation_Type, Position) %>%
    summarize(
      week       = week[which.min(abs(week - median(week)))],
      percentage = percentage[which.min(abs(week - median(week)))],
      Product    = first(Product),
      .groups    = "drop"
    )
  
  pal <- c(
    "#1f77b4", "#d99255", "#448744", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#7d485f", "#c7c7c7", "#8a7322", "#2a7785"
  )

  # Plot: time-series with one line per Position, labels spread via repel with one line per Position, labels spread via repel
  p <- ggplot(df, aes(x = week, y = percentage, color = Position)) +
    geom_line(aes(group = Position), alpha = 0.5) +
    geom_point(size = 1) +
    geom_text_repel(
      data          = label_df,
      aes(label = str_wrap(Product, width = 10)),
      size          = 2,
      box.padding   = 0.1,
      point.padding = 0.1,
      segment.size  = 0.1
    ) +
    facet_wrap(~ Mutation_Type) +
    scale_color_manual(values = pal, name = "Position") +
    labs(
      title    = paste0(samp, " — ", sp),
      subtitle = paste0("Top ", length(top_positions), "/",
                        n_distinct(df$Position),
                        " by p-value (<= ", signif(cutoff_p, 2), ")"),
      x = "Week",
      y = "Mutation Frequency"
    ) +
    custom_theme
  
  # Save plot with dynamic filename
  fname <- sprintf("%s_%s_top%d_p%s_time_series.png",
                   samp, sp,
                   length(top_positions),
                   gsub("0\\.", "",
                        format(signif(cutoff_p, 2), scientific = FALSE)))
  ggsave(
    file.path(fig_dir, fname),
    plot   = p,
    height = 6, width = 12, units = "in", dpi = 300
  )
}

# -----------------------------
# End of iterative analysis
# -----------------------------

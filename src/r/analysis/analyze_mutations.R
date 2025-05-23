# analyze_mutations_refactored.R -------------------------------------------------
# Date: 2025‑05‑22
#
# This refactor consumes the **Parquet summaries** emitted by
# `summarise_mutations.py` (v5) and produces publication‑ready statistics
# and figures.  It is fully tidyverse‑centric, arrow‑based (no readr bottlenecks),
# and broken into modular functions so you can source() pieces in an RMarkdown.
# ------------------------------------------------------------------------------

# ── Libraries ────────────────────────────────────────────────────────────────
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  arrow,            # efficient Parquet loading
  dplyr, tidyr,     # data wrangling
  ggplot2, ggrepel, # plotting
  ggridges,         # ridge plots for trajectories
  patchwork,        # plot assembly
  janitor,          # clean_names()
  lubridate,        # date helpers (optional)
  stringr,
  duckdb,
  purrr
)

#family = "Arial", 
theme_pub <- function() {
  theme_classic() +
  theme(
    text = element_text(family = "Arial", size = 12, color = "black"),  # Customize font, size, and color
    axis.title.x = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 14),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_line(linewidth = 0.75, color = "black"),
    axis.ticks.length = unit(0.3, "cm"),  # Make tick marks longer
    axis.line = element_line(linewidth = 0.75, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 12, color = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    legend.position = "none",
    strip.background = element_blank(),  # Remove the box from facet wrap titles
    strip.text = element_text(size = 12, color = "black")  # Customize facet title text
  )
}

# ── Helper: CLI args ---------------------------------------------------------
cli_opts <- list(
  in_dir      = "/Users/reneeoles/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Research/Strain_Evolution/iHMP/output",          # where the Parquet files live
  meta_file   = "../metadata/hmp2_metadata_2018-08-20.csv",
  out_dir     = "figures",            # where to save PDFs/PNGs
  dpi         = 320
)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  for (a in args) eval(parse(text = a))  # allow --in_dir="…" style
}
cli_opts <- modifyList(cli_opts, mget(names(cli_opts)))
if (!dir.exists(cli_opts$out_dir)) dir.create(cli_opts$out_dir, recursive = TRUE)

# ── Load metadata ------------------------------------------------------------
meta <- readr::read_csv(cli_opts$meta_file) %>%
  rename(patient_id = `Participant ID`) %>%
  select(patient_id, diagnosis, week_num, hbi, sccai, "fecalcal")

# ── Arrow datasets -----------------------------------------------------------
open_pq <- function(name) arrow::open_dataset(file.path(cli_opts$in_dir, paste0(name, ".parquet")))

gene_traj  <- open_pq("gene_traj")
bin_dnds   <- open_pq("bin_dnds")
patient_wk <- open_pq("patient_week")
snp_stats  <- open_pq("snp_stats")

# ── 1. Allele‑frequency trajectories by diagnosis ---------------------------
plot_trajectories <- function() {
  traj_df <- arrow::to_duckdb(snp_stats) %>%
    filter(p_value < 0.05) %>%
    select(patient_id, chromosome, position, slope, p_value) %>%
    collect() %>%
    left_join(meta, by = "patient_id")

  ggplot(traj_df, aes(x = slope, y = diagnosis, fill = diagnosis)) +
    ggridges::geom_density_ridges(alpha = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Distribution of per SNV slopes", x = "Linear trajectory slope", y = NULL) +
    theme_pub()
}

# ── 2. Gene volcano plot -----------------------------------------------------
plot_gene_volcano <- function() {
  volc <- gene_traj %>%
    group_by(gene) %>%
    summarise(min_p = min(min_p), mean_abs_slope = mean(mean_abs_slope, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(log_p = -log10(min_p)) %>%
    collect()

  ggplot(volc, aes(mean_abs_slope, log_p)) +
    geom_point(alpha = 0.6) +
    ggrepel::geom_text_repel(data = subset(volc, log_p > 10 | mean_abs_slope > 0.3), aes(label = gene), size = 3) +
    labs(title = "Gene‑level signal volcano", x = "Mean |slope|", y = "‑log10(min P)") +
    theme_pub()
}

# ── 3. Bin dN/dS summary -----------------------------------------------------
plot_bin_dnds <- function() {
  dnds <- bin_dnds %>%
    mutate(log_dnds = log10(dnds+0.001)) %>%
    filter(is.finite(log_dnds)) %>%
    collect() %>%
    left_join(meta, by = "patient_id")

  ggplot(dnds[dnds$log_dnds>-3,], aes(log_dnds)) +
    geom_histogram(bins = 60, fill = "steelblue", colour = "black", alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = "Bin‑level dN/dS (log10)", x = "log10(dN/dS)", y = "Bins") +
    facet_wrap(~diagnosis) +
    theme_pub()

  ggplot(dnds[dnds$log_dnds>-3,], aes(y=log_dnds,x=diagnosis)) +
    geom_boxplot() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = "Bin‑level dN/dS (log10)", x = "log10(dN/dS)", y = "Bins") +
    theme_pub()
}

# ── 4. Mutation burden over time --------------------------------------------
plot_burden <- function() {
  burden <- patient_wk %>%
    mutate(week_num = as.numeric(week_num)) %>%  # coerce to numeric
    collect() %>%
    left_join(meta %>% mutate(week_num = as.numeric(week_num)), by = c("patient_id", "week_num"))
  
  ggplot(burden, aes(x = week_num, y = n_snvs, group = patient_id, colour = diagnosis)) +
    geom_line(alpha = 0.5) +
    geom_smooth(method = "loess", se = FALSE, size = 1.2) +
    scale_colour_brewer(palette = "Dark2") +
    labs(title = "SNV burden per patient over time", x = "Week", y = "Number of SNVs") +
    theme_pub()
}

# ── Save all plots -----------------------------------------------------------
save_plot <- function(p, name, w = 6, h = 4) {
  ggsave(file.path(cli_opts$out_dir, paste0(name, ".pdf")),
         plot = p, width = w, height = h, device = "pdf")
  ggsave(file.path(cli_opts$out_dir, paste0(name, ".png")),
         plot = p, width = w, height = h, dpi = cli_opts$dpi)
}

message("📊 Generating plots …")
plots <- list(
  traj      = plot_trajectories(),
  volcano   = plot_gene_volcano(),
  bin_dnds  = plot_bin_dnds(),
  burden    = plot_burden()
)

walk2(plots, names(plots), save_plot)

message("✅ All figures written to ", cli_opts$out_dir)

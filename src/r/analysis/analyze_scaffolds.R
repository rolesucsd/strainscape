# analyze_scaffolds_refactored.R
# ---------------------------------
# Author: Strainscape Assistant (ChatGPT)
# Last updated: 2025‑05‑22
#
# Purpose
# -------
# Provide a fully‑annotated, publication‑ready analysis pipeline for combined
# InStrain scaffold summaries that 
#   * merges sample‑level scaffold data with iHMP metadata,
#   * produces longitudinal summaries per patient and per week,
#   * reports bin‑level summaries (while keeping bins sample‑specific), and
#   * generates high‑quality graphical outputs suitable for a manuscript.
#
# Usage (minimal example)
# ----------------------
# source("analyze_scaffolds_refactored.R")
# results <- analyze_scaffolds(
#              scaffold_file  = "combined_processed_scaffolds_sample.txt",
#              metadata_file  = "hmp2_metadata_2018-08-20.csv",
#              out_dir        = "output"  # will be created if it doesn't exist
#            )
#
# R/CRAN packages required
# -----------------------
# tidyverse  >= 2.0.0
# janitor    >= 2.2.0   (clean names)
# patchwork  >= 1.2.0   (plot assembly)
# ggpubr     >= 0.6.0   (non‑parametric stats on plots)
# scales     >= 1.3.0
# ggrepel    >= 0.9.5   (optional – label repelling)
#
# ---------------------------------------------------------------------------
# 1. Helper functions --------------------------------------------------------
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(patchwork)
  library(ggpubr)
  library(scales)
})

#' Read scaffold table (TSV) and clean column names.
#' @param path  Character. Path to combined_processed_scaffolds file.
#' @return tibble
read_scaffolds <- function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    clean_names() %>%
    mutate(sample = as.character(sample))
}

#' Read HMP2 metadata (CSV) and subset key columns.
#' @param path Character. Metadata path.
#' @param keep Default c("External.ID","diagnosis","Participant.ID","week.num")
#' @return tibble with renamed keys (sample, patient_id, diagnosis, week)
read_metadata <- function(path,
                          keep = c("External.ID", "diagnosis", "Participant ID", "week_num")) {
  read_csv(path, show_col_types = FALSE) %>%
    select(all_of(keep)) %>%
    janitor::clean_names() %>%
    rename(sample = external_id,
           patient_id = participant_id,
           week = week_num)
}

#' Merge scaffold & metadata tables.
#' @inheritParams read_scaffolds
#' @inheritParams read_metadata
#' @return joined tibble
merge_scaffold_metadata <- function(scaffolds, metadata) {
  left_join(scaffolds, metadata, by = "sample")
}

#' Longitudinal summaries per patient x week.
#' Generates three summary tables (coverage, breadth, diversity) each in wide
#' and long forms for flexible downstream plotting.
patient_week_summaries <- function(data) {
  long <- data %>%
    group_by(patient_id, week) %>%
    summarise(n_scaffolds   = n(),
              mean_cov      = mean(coverage, na.rm = TRUE),
              sd_cov        = sd(coverage,   na.rm = TRUE),
              mean_breadth  = mean(breadth,   na.rm = TRUE),
              mean_div      = mean(nucl_diversity, na.rm = TRUE),
              .groups = "drop")

  wide <- long %>%
    pivot_wider(names_from = week,
                values_from = c(n_scaffolds, mean_cov, mean_breadth, mean_div))

  list(long = long, wide = wide)
}

#' Bin‑level summaries per patient (bins are sample‑specific!)
#' We therefore summarise within each patient *and* sample.
#' @return tibble of mean metrics per (patient, week, bin), including completeness and contamination
patient_bin_summaries <- function(data) {
  data %>%
    group_by(patient_id, week, bin) %>%
    summarise(mean_cov   = mean(coverage, na.rm = TRUE),
              mean_bread = mean(breadth,  na.rm = TRUE),
              mean_div   = mean(nucl_diversity, na.rm = TRUE),
              completeness = first(Completeness),
              contamination = first(Contamination),
              .groups = "drop")
}

#' Save a tidy CSV to out_dir with a descriptive filename.
write_out <- function(df, out_dir, fname) {
  write_csv(df, file.path(out_dir, fname))
}

# ---------------------------------------------------------------------------
# 2. Plotting utilities ------------------------------------------------------
# ---------------------------------------------------------------------------

#' Generic ggplot theme for publication.
plot_theme <- function() {
  theme_classic() +
  theme(
    text = element_text(family = "sans", size = 12, color = "black"),  # Changed from Arial to sans
    axis.title.x = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 14),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_line(size = 0.75, color = "black"),
    axis.ticks.length = unit(0.3, "cm"),  # Make tick marks longer
    axis.line = element_line(size = 0.75, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 12, color = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    legend.position = "none",
    strip.background = element_blank(),  # Remove the box from facet wrap titles
    strip.text = element_text(size = 12, color = "black")  # Customize facet title text
  )
}

#' Line plot of mean coverage over weeks for each patient.
plot_cov_trend <- function(df) {
  ggplot(df, aes(week, mean_cov, group = patient_id, colour = patient_id)) +
    geom_line() + geom_point(size = 2) +
    scale_y_continuous("Mean scaffold coverage", trans = "log10",
                       labels = label_number()) +
    scale_x_continuous("Study week", breaks = pretty_breaks()) +
    labs(title = "Longitudinal coverage per patient") +
    plot_theme()
}

#' Violin + jitter plot of bin diversity by diagnosis.
plot_bin_div_by_dx <- function(bin_stats) {
  ggplot(bin_stats, aes(diagnosis, log10(mean_div), fill = diagnosis)) +
    geom_boxplot(trim = FALSE) +
  geom_hline(yintercept = c(-3,-2), linetype = "dashed", color = "grey") +
  stat_compare_means(comparisons = list(c("nonIBD","UC"),
                                        c("nonIBD","CD"),
                                        c("UC","CD")),
                     method = "wilcox.test", label = "p.format") +
  stat_compare_means(method = "kruskal.test", label.y = -0.5) +  # Add Kruskal-Wallis test
  labs(y = "Log10 mean nucleotide diversity",
         title = "Bin level diversity by diagnosis") +
  scale_fill_manual(values = c("UC"="#bb9df5",
                               "nonIBD"="#71d987",
                               "CD"="#ffc169")) +
  plot_theme()
}

# ---------------------------------------------------------------------------
# 3. Main user‑facing wrapper ------------------------------------------------
# ---------------------------------------------------------------------------

#' Analyze iHMP scaffold output with metadata context.
#'
#' @param scaffold_file Path to combined_processed_scaffolds file (tsv).
#' @param metadata_file Path to hmp2 metadata CSV.
#' @param out_dir       Directory for all outputs (created if missing).
#' @return named list of tibbles + ggplot objects.
#' @export
analyze_scaffolds <- function(scaffold_file, metadata_file, out_dir = "output") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # ---------- I/O ----------
  scaffolds <- read_scaffolds(scaffold_file)
  meta      <- read_metadata(metadata_file)
  dat       <- merge_scaffold_metadata(scaffolds, meta)

  # ---------- Summaries ----------
  pw  <- patient_week_summaries(dat)
  pb  <- patient_bin_summaries(dat)

  # ---------- Write summaries ----------
  write_out(pw$long, out_dir, "patient_week_long.csv")
  write_out(pw$wide, out_dir, "patient_week_wide.csv")
  write_out(pb,       out_dir, "patient_bin_summary.csv")

  # ---------- Figures ----------
  cov_plot <- plot_cov_trend(pw$long)
  ggsave(file.path(out_dir, "coverage_trend.png"),  cov_plot, width = 7, height = 4)

  bin_plot <- plot_bin_div_by_dx(pb %>%
                                   left_join(meta, by = "patient_id"))
  ggsave(file.path(out_dir, "bin_diversity_diagnosis.png"), bin_plot, width = 6, height = 6)

  # ---------- Return ----------
  list(raw            = dat,
       patient_week   = pw,
       patient_bin    = pb,
       figures        = list(coverage_trend = cov_plot,
                             bin_diversity  = bin_plot))
}

# ---------------------------------------------------------------------------
# 4. Command‑line interface --------------------------------------------------
# ---------------------------------------------------------------------------
if (identical(Sys.getenv("R_SCRIPT_RUN"), "TRUE")) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Usage: R_SCRIPT_RUN=TRUE Rscript analyze_scaffolds_refactored.R <scaffold.tsv> <metadata.csv> [out_dir]",
         call. = FALSE)
  }
  scaffold_file <- args[[1]]
  metadata_file <- args[[2]]
  out_dir       <- ifelse(length(args) >= 3, args[[3]], "output")
  analyze_scaffolds(scaffold_file, metadata_file, out_dir)
}

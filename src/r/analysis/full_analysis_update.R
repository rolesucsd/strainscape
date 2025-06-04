# ── packages -----------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  arrow, dplyr, tidyr, tibble, purrr, janitor,
  DESeq2, broom, stringr,
  ggplot2, ggrepel, patchwork, scales, ggpub
)

# ── ggplot theme -------------------------------------------------------------
theme_pub  <- function() {
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
theme_set(theme_pub())

# ── paths ---------------------------------------------------
DIR  <- "/Users/reneeoles/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Research/Strain_Evolution/iHMP/output"
META <- "../metadata/hmp2_metadata_2018-08-20.csv"
BINS <- "patient_bin_summary.csv"
OUT  <- "figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ════════════════════════════════════════════════════════════════════════════
# 1.  Metadata
# ════════════════════════════════════════════════════════════════════════════
meta <- read.csv(META) |>
  janitor::clean_names() |>
  dplyr::rename(patient_id = participant_id) |>
  dplyr::mutate(group = case_when(
    diagnosis == "nonIBD" ~ "nonIBD",
    diagnosis == "UC"     ~ "UC",
    diagnosis == "CD"     ~ "CD",
    TRUE                  ~ NA_character_
  )) |>
  dplyr::select(patient_id, group) |>
  dplyr::distinct()

# ════════════════════════════════════════════════════════════════════════════
# 2.  Nucleotide diversity
# ════════════════════════════════════════════════════════════════════════════
bins_summary <- read_feather(file.path(DIR, "prep/bins_summary.feather"))
bins_summary <- bins_summary %>%
  mutate(group = factor(group, levels = c("nonIBD", "UC", "CD")))
ggplot(bins_summary, aes(group, log10(mean_div), fill = group)) +
  geom_boxplot(trim = FALSE) +
  geom_hline(yintercept = c(-3,-2), linetype = "dashed", color = "grey") +
  stat_compare_means(comparisons = list(c("nonIBD","UC"),
                                        c("nonIBD","CD"),
                                        c("UC","CD")),
                     method = "wilcox.test", label = "p.format") +
  labs(y = "Log10 mean nucleotide diversity",
       title = "Bin level diversity by diagnosis") +
  scale_fill_manual(values = c("UC"="#bb9df5",
                               "nonIBD"="#71d987",
                               "CD"="#ffc169"))
ggsave(file.path(OUT, "bin_nuc_div.png"), width = 2.5, height = 4)

## ─────────────────────────────────────────────────────────────────────────
## 3.  SNVs per bin
## ─────────────────────────────────────────────────────────────────────────
snvs_per_bin <- read_feather(file.path(DIR, "prep/snvs_per_bin.feather"))
ggplot(snvs_per_bin, aes(x=`TRUE`)) +
  geom_histogram(bins = 8, color="white") +
  scale_x_log10() +
  labs(x = "Log10 No. of Sweeps", y = "No. of isolates")
ggsave(file.path(OUT, "sweep_histogram.png"), width = 2.5, height = 4)

## ─────────────────────────────────────────────────────────────────────────
## 4.  Sweeps per bin
## ─────────────────────────────────────────────────────────────────────────
bins_summary_filter <- read_feather(file.path(DIR, "prep/bins_summary_filter.feather"))
# filter out the scaffold/bins output to only include isolates with less than/equal to 1000 snvs
table(bins_summary_filter$group)
#nonIBD     UC     CD 
#659    565    947 
ggplot(bins_summary_filter, aes(y=-log10(prop), x=group, fill=group))+
  geom_boxplot()+
  stat_compare_means(comparisons = list(c("nonIBD","UC"),
                                        c("nonIBD","CD"),
                                        c("UC","CD")),
                     method = "wilcox.test", label = "p.format") +
  #  stat_compare_means(method = "kruskal.test", label.y = -0.5) +  # Add Kruskal-Wallis test
  labs(y = "Log10 sweeping mutations/total", x="") +
  scale_fill_manual(values = c("UC"="#bb9df5",
                               "nonIBD"="#71d987",
                               "CD"="#ffc169"))
ggsave(file.path(OUT, "sweep_proportion.png"), width = 2.5, height = 4)


## ─────────────────────────────────────────────────────────────────────────
## 5.  Sweep enrichment analyses
## ─────────────────────────────────────────────────────────────────────────
enrichment <- read_feather(file.path(DIR, "prep/enrichment_results.feather"))

# Now you can graph each individually
df <- enrichment %>%
  filter(subset %like% "1k")

df <- df %>%
  mutate(sig  = enrichment > 2 & p_value < 0.05)
    
ggplot(df, aes(logE, logP, colour = sig)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept =  log10(2),     linetype = "dashed") +
  geom_point(size = 3, alpha = 0.7) +
  ggrepel::geom_text_repel(data = dplyr::filter(df, sig),aes(label = feature), size = 3,max.overlaps = 100) +
  scale_colour_manual(values = c("TRUE" = "#ffc169", "FALSE" = "grey70")) +
  labs(
    x     = bquote(log[10]~"enrichment"),
    y     = expression(-log[10]~p-value)
  ) +
  xlim(0, NA) +
  theme_pub()

fname <- file.path(
  OUT,
  paste0("volcano_", "snv_sweeps_ec_terms_CD_vs_nonIBD", ".png")
)
ggsave(fname, width = 3.5, height = 3.5)
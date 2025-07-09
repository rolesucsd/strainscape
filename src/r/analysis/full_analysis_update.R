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

bins_summary <- bins_summary %>%
  mutate(group2 = if_else(group == "nonIBD", "nonIBD", "IBD"))
tt <- t.test(mean_div ~ group2, data = bins_summary)
tidy(tt)

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
  filter(subset %in% "snv_sweeps") %>%
  filter(feature_type %in% "product")


df <- read.delim("/Users/reneeoles/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Research/Strain_Evolution/iHMP/output/snv_sweeps_genes_vs_nonIBD.txt")
df <- df %>%
  mutate(sig = case_when(
    enrichment > 2 & FDR < 0.1 & grepl("UC", comparison) ~ "UC",
    enrichment > 2 & FDR < 0.1 & grepl("CD", comparison) ~ "CD",
    TRUE ~ "nonsig"
  ))
ggplot(df, aes(x = log10(enrichment), y = -log10(FDR))) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
  geom_vline(xintercept =  log10(2),     linetype = "dashed") +
  geom_point(
    data = df[df$sig == "nonsig",],
    aes(colour = sig),
    size = 1, alpha = 0.7
  ) +
  geom_point(
    data = subset(df, sig != "nonsig"),
    aes(colour = sig),
    size = 2, alpha = 0.7
  ) +
  geom_text_repel(
    data = df[df$sig != "nonsig",],
    aes(label       = Annotation),
    size           = 6,
    color          = "black",
    max.overlaps   = 50
  ) +
  scale_color_manual(values = c(
    "UC"     = "#877b9e",
    "nonsig" = "grey",
    "CD"     = "#9c7d52"
  )) +
  labs(
    x = bquote(log[10]~"enrichment"),
    y = expression(-log[10]~"p-value")
  ) +
  xlim(0, NA) +
  theme_pub()
fname <- file.path(
  OUT,
  paste0("volcano_", "snv_sweeps_gene_terms_vs_nonIBD", ".png")
)
ggsave(fname, width = 6, height = 5.5)


## ─────────────────────────────────────────────────────────────────────────
## 6.  CPX
## ─────────────────────────────────────────────────────────────────────────
meta_fecalcal <- read.csv("/Users/reneeoles/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Research/Strain_Evolution/iHMP/metadata/hmp2_metadata_2018-08-20.csv")
colnames(meta_fecalcal)[3] <- "patient_id"
meta_fecalcal <- meta_fecalcal %>%  
  select(patient_id, week_num, fecalcal) %>%
  filter(complete.cases(.))

snv_raw <- read_tsv("snv_vs_cpx_pearson.tsv", col_types = cols())
num_cols         <- names(snv_raw) %>% keep(~ grepl("^\\d+$", .x))

snv_meta <- snv_raw %>% 
  left_join(select(bins_summary_filter, patient_id, bin, 9:11), by = c("patient_id","bin")) %>% 
  left_join(select(bin_indentify, patient_id, bin, classification)) %>%
  filter(abs(r) > 0.99)
write.table(snv_meta, "snv_vs_cpx_pearson_r099.txt", row.names = F, sep="\t", quote = F)

snv_numeric      <- snv_meta %>% select(patient_id, bin, position, gene, all_of(num_cols))
snv_meta         <- snv_meta %>% 
  select(-all_of(num_cols)) %>%
  separate_rows(gene, sep = "/") %>%
  filter(!is.na(gene)) %>%
  filter(!(gene %in% "nan"))
  snv_meta_select <- snv_meta %>%
  filter(gene %in% c(oxidative_genes, "fecR", "hemW", "fixB", "cydB", "hydF", "lpdA", "adhE"))

snv_long <- snv_numeric %>%
  filter(gene %in% c(oxidative_genes, "fecR", "hemW", "fixB", "cydB", "hydF", "lpdA", "adhE")) %>%
#  filter(patient_id == pid_target,
#         position    == pos_target,
#         bin  == chr_target) %>%
  pivot_longer(
    cols = where(is.numeric) & !c(patient_id, bin, position),
    names_to  = "week_num",
    values_to = "Freq",
    values_drop_na = TRUE
  ) %>%
  mutate(week_num = as.integer(week_num))

snv_cpx_with_cpx <- snv_long %>%
  left_join(meta_fecalcal, by = c("patient_id", "week_num"))

snv_cpx_with_cpx %>%
  ggplot(aes(y=Freq, x=fecalcal, color = gene)) +
    geom_point()+
    geom_line() +
    theme_pub() +
  theme(legend.position = "bottom")
fname <- file.path(
  OUT,
  paste0("snv_ox_genes_fecalcal", ".png")
)
ggsave(fname, width = 4, height = 4)


## ─────────────────────────────────────────────────────────────────────────
## Look at rates of substitutions among sweep events
## ─────────────────────────────────────────────────────────────────────────
sweeps <- sweeps %>%
  mutate(substitution = paste0(ref_base, ">", new_base)) %>%
  left_join(meta) %>%
  left_join(bin_indentify) %>%
  mutate(species = sub(".*s__", "", classification)) %>%
  mutate(
    GeneType = if_else(gene %in% c(oxidative_genes), "ros", "other")
  )

sweep_counts <- sweeps %>%
  dplyr::count(group, species, patient_id, bin, substitutions, name = "count")

all_subs      <- sort(unique(sweep_counts$substitutions))

sweeps_prop <- sweep_counts %>% 
  group_by(group, species, patient_id, bin) %>%
  complete(
    substitutions  = all_subs,
    fill          = list(count = 0)
  ) %>% 
  ungroup() %>%
  group_by(group, species, patient_id, bin) %>%
  mutate(prop = count / sum(count)) %>% 
  ungroup()

# One plot not differentiating ROS or species
sweeps_prop %>%
  mutate(group = factor(group, levels = c("nonIBD", "UC", "CD"))) %>%
  ggplot(aes(x=group, y=prop)) +
  stat_summary(
    fun    = mean,
    geom   = "col",
    position = position_dodge(width = 0.8),
    width    = 0.7,
    alpha    = 0.8,
    fill = "grey70"
  ) +  
  scale_y_continuous(labels = percent_format(1)) +
  labs(
    x = "Group",
    y = "Percentage of mutations"
  ) +
  facet_wrap(~substitutions, scales="free_y") +
  custom_theme  +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
fname <- file.path(
  OUT,
  paste0("sweeps_substitution_rates", ".png")
)
ggsave(fname, width = 5, height = 5)

sweeps_prop %>%
  filter(substitutions %like% "C>A") %>%
  group_by(species) %>%
  filter(n() >= 20) %>%
  ungroup()%>%
  mutate(group = factor(group, levels = c("nonIBD", "UC", "CD"))) %>%
  ggplot(aes(
    x     = group,
    y     = prop,
    fill = group
  )) +
  stat_summary(
    fun    = mean,
    geom   = "col",
    position = position_dodge(width = 0.8),
    width    = 0.7,
    alpha    = 0.8
  ) +  
  facet_wrap(~ species) +
  scale_y_continuous(labels = percent_format(1)) +
  labs(
    x = "Group",
    y = "Proportion of C→A/A→C mutations",
    fill = "Gene type"
  ) +
  scale_fill_manual(values = c("UC"="#bb9df5",
                               "nonIBD"="#71d987",
                               "CD"="#ffc169")) +
  custom_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) 
fname <- file.path(
    OUT,
    paste0("sweeps_substitution_rates_species", ".png")
  )
ggsave(fname, width = 10, height = 9)


sweep_counts <- sweeps %>%
  mutate(
    GeneType = if_else(gene %in% c(oxidative_genes), "ros", "other")
  ) %>%
  dplyr::count(group, GeneType, species, patient_id, bin, substitutions, name = "count")

all_subs      <- sort(unique(sweep_counts$substitutions))

sweeps_prop <- sweep_counts %>% 
  group_by(group, species, patient_id, bin, GeneType) %>%
  complete(
    substitutions  = all_subs,
    fill          = list(count = 0)
  ) %>% 
  ungroup() %>%
  group_by(group, species, patient_id, bin, GeneType) %>%
  mutate(prop = count / sum(count)) %>% 
  ungroup() 


# One plot not differentiating ROS or species
sweeps_prop %>%
  mutate(group = factor(group, levels = c("nonIBD", "UC", "CD"))) %>%
  ggplot(aes(x=group, y=prop, fill=GeneType)) +
  stat_summary(
    fun    = mean,
    geom   = "col",
    position = position_dodge(width = 0.8),
    width    = 0.7,
    alpha    = 0.8
  ) +  
  scale_fill_manual(values = c(ros   = "blue",
                                other = "grey70")) +
  scale_y_continuous(labels = percent_format(1)) +
  labs(
    x = "Group",
    y = "Percentage of mutations"
  ) +
  facet_wrap(~substitutions, scales="free_y") +
  custom_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) 
fname <- file.path(
  OUT,
  paste0("sweeps_substitution_rates_ros", ".png")
)
ggsave(fname, width = 5, height = 5)


sweeps_prop %>%
  filter(substitutions %like% "C>A") %>%
  group_by(species) %>%
  filter(n() >= 10) %>%
  ungroup()%>%
  ggplot(aes(
    x     = group,
    y     = prop,
    fill  = GeneType,
  )) +
  stat_summary(
    fun    = mean,
    geom   = "col",
    position = position_dodge(width = 0.8),
    width    = 0.7,
    alpha    = 0.8
  ) +  
#  geom_jitter(size = 1.5, alpha = 0.7, position = position_jitterdodge(dodge.width  = 0.8,jitter.width = 0.05)) +
  facet_wrap(~ species, scales = "free_y") +
  scale_y_continuous(labels = percent_format(1)) +
  scale_fill_manual(values = c(ros   = "blue",
                               other = "grey70")) +
  labs(
    x = "Group",
    y = "Proportion of C→A/A→C mutations",
    fill = "Gene type"
  ) +
  custom_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


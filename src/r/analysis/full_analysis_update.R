# ── packages -----------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  arrow, dplyr, tidyr, tibble, purrr, janitor,
  DESeq2, broom, stringr,
  ggplot2, ggrepel, patchwork, scales, ggpub
)

# ── ggplot theme -------------------------------------------------------------
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

bins <- read.csv(BINS)
bins_summary <- bins %>%
  group_by(patient_id, bin) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
    .groups = "drop"
  )
bins_summary <- left_join(bins_summary, meta)

bins_summary <- bins_summary %>%
  mutate(group = factor(group, levels = c("nonIBD", "UC", "CD")))
library(ggpubr)
ggplot(bins_summary, aes(group, log10(mean_div), fill = group)) +
  geom_boxplot(trim = FALSE) +
  geom_hline(yintercept = c(-3,-2), linetype = "dashed", color = "grey") +
  stat_compare_means(comparisons = list(c("nonIBD","UC"),
                                        c("nonIBD","CD"),
                                        c("UC","CD")),
                     method = "wilcox.test", label = "p.format") +
#  stat_compare_means(method = "kruskal.test", label.y = -0.5) +  # Add Kruskal-Wallis test
  labs(y = "Log10 mean nucleotide diversity",
       title = "Bin level diversity by diagnosis") +
  scale_fill_manual(values = c("UC"="#bb9df5",
                               "nonIBD"="#71d987",
                               "CD"="#ffc169"))
ggsave(file.path(OUT, "bin_nuc_div.png"), width = 2.5, height = 4)


## ─────────────────────────────────────────────────────────────────────────
## 2.  Load SNVs
## ─────────────────────────────────────────────────────────────────────────
snvs <- arrow::open_dataset(file.path(DIR, "all_snvs")) %>%
  collect()

snvs <- snvs %>%
  mutate(is_sweep = abs(min_freq-max_freq) >= 0.7)

# look at just the snvs with selective sweep 
sweeps <- snvs |> filter(is_sweep == TRUE)

# Calculate the number of snvs per isolate (patiebnt + bin combination)
snvs_per_bin <- snvs %>%
  group_by(patient_id, bin, is_sweep) %>%
  summarise(count = n(), .groups = "drop")
snvs_per_bin <- snvs_per_bin %>%
  pivot_wider(
    names_from = is_sweep,
    values_from = count,
    values_fill = 0
  )

ggplot(snvs_per_bin, aes(x=`TRUE`)) +
  geom_histogram(bins = 8, color="white") +
  scale_x_log10() +
  labs(x = "Log10 No. of Sweeps", y = "No. of isolates")
ggsave(file.path(OUT, "sweep_histogram.png"), width = 2.5, height = 4)

bins_summary_filter <- inner_join(bins_summary, snvs_per_bin)
bins_summary_filter <- bins_summary_filter[bins_summary_filter$`TRUE` <= 1000,]
bins_summary_filter <- na.omit(bins_summary_filter)
bins_summary_filter$prop <- bins_summary_filter$`TRUE`/bins_summary_filter$`FALSE`
bins_summary_filter$total <- bins_summary_filter$`TRUE` + bins_summary_filter$`FALSE`
# filter out the scaffold/bins output to only include isolates with less than/equal to 1000 snvs
table(bins_summary_filter$group)
#nonIBD     UC     CD 
#660    569    956 

library(ggpubr)
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


# filter the snvs output (goes from 9 million to 279349)
sweeps_filtered <- left_join(bins_summary_filter, sweeps)

first3  <- names(sweeps_filtered)[c(1:10,31:48,88:92)]
num_cols <- names(sweeps_filtered)[
  grepl("^\\d+$", names(sweeps_filtered))]
num_cols <- num_cols[order(as.numeric(num_cols))]
others   <- setdiff(names(sweeps_filtered), c(first3, num_cols))
sweeps_filtered <- sweeps_filtered[c(first3, num_cols, others)]

#write.table(sweeps_filtered, "sweeps_filtered.txt", row.names = F, sep="\t", quote=F)

# Input iterations
# 1. snvs - all sweeps, no filter
snv_sweeps <- sweeps %>%
  left_join(meta) %>%
  distinct(patient_id, bin, group, product, gene, kegg_terms, go_terms, ec_terms, dbxrefs, locus_tag)
# 2. snvs - all sweeps, no filte
snv_sweeps_nosilent <- sweeps[!(sweeps$mutation_type %in% "Silent"),] %>%
  left_join(meta) %>%
  distinct(patient_id, bin, group, product, gene, kegg_terms, go_terms, ec_terms, dbxrefs, locus_tag)
# 3. snvs - filter by sweeps first, then filter by bins with greater than 1k snvs 
snv_sweeps_1k <- sweeps_filtered %>%
  distinct(patient_id, bin, group, product, gene, kegg_terms, go_terms, ec_terms, dbxrefs, locus_tag)
# 4. snvs - remove silent mutations 
snv_sweeps_1k_nosilent <- sweeps_filtered[!(sweeps_filtered$mutation_type %in% "Silent"),] %>%
  distinct(patient_id, bin, group, product, gene, kegg_terms, go_terms, ec_terms, dbxrefs, locus_tag)

## ─────────────────────────────────────────────────────────────────────────
## 3.  Now I want to run this section both with the column Gene and Product for each of the 4 input
## ─────────────────────────────────────────────────────────────────────────
fisher_enrich <- function(df, col) {
  expanded <- df %>%
    mutate(row_id = row_number()) %>%
    group_by(row_id) %>%
    group_modify(~ {
      col_sym <- sym(col)
      value <- .x[[col]]
      
      # Check for intergenic and `/` present
      if (.x$mutation_type == "intergenic" && str_detect(value, "/")) {
        parts <- unlist(str_split(value, "/"))
        n <- length(parts)
        if (n < 2) return(.x)
        
        mid <- ceiling(n / 2)
        new1 <- paste(parts[1:mid], collapse = "/")
        new2 <- paste(parts[(mid + 1):n], collapse = "/")
        
        .x1 <- .x
        .x2 <- .x
        .x1[[col]] <- new1
        .x2[[col]] <- new2
        bind_rows(.x1, .x2)
      } else {
        .x
      }
    }) %>%
    ungroup() %>%
    {
      if (str_detect(col, "term")) {
        tidyr::separate_rows(., !!sym(col), sep = ",", convert = FALSE)
      } else {
        .
      }
    } %>%
    mutate(!!sym(col) := str_trim(.data[[col]])) %>%
    filter(!is.na(.data[[col]]),
           .data[[col]] != "nan",
           .data[[col]] != "hypothetical protein") %>%
    select(-row_id)
  
  # STEP C – per-feature, per-group bin counts
  feature_counts <- expanded %>%
    group_by(!!sym(col), group) %>%
    summarise(count = n(), .groups = "drop")
  
  # STEP D – group sizes (# of (patient,bin) per diagnosis group)
  group_sizes <- expanded %>%
    dplyr::distinct(patient_id, bin, group) %>%
    dplyr::count(group, name = "total_bins") %>%
    tibble::deframe()
  
    # STEP E – Fisher for every pair of groups  (unchanged)
  pairs <- list(
    CD_vs_nonIBD = c("CD",      "nonIBD"),
    CD_vs_UC     = c("CD",      "UC"),
    UC_vs_nonIBD = c("UC",      "nonIBD")
  )
  
  purrr::imap_dfr(pairs, function(groups, cname) {
    g1 <- groups[1]; g2 <- groups[2]
    
    feature_counts %>%
      filter(group %in% groups) %>%
      pivot_wider(names_from = group, values_from = count, values_fill = 0) %>%
      mutate(
        total1 = group_sizes[[g1]],
        total2 = group_sizes[[g2]],
        n1     = .[[g1]],
        n2     = .[[g2]]
      ) %>%
      rowwise() %>%
      mutate(
        ok  = all(is.finite(c(n1, n2, total1, total2))) &&
          all(c(n1, total1 - n1, n2, total2 - n2) >= 0),
        ft  = if (ok) list(fisher.test(matrix(c(n1,
                                                total1 - n1,
                                                n2,
                                                total2 - n2), 2))) else list(NULL),
        p_value    = if (!is.null(ft)) ft$p.value    else NA_real_,
        odds_ratio = if (!is.null(ft)) unname(ft$estimate) else NA_real_,
        enrichment = ((n1 + 0.5) / total1) / ((n2 + 0.5) / total2)
      ) %>%
      ungroup() %>%
      select(-ft, -ok) %>%
      transmute(
        comparison     = cname,
        feature        = .data[[col]],
        !!paste0("n_",     g1) := n1,
        !!paste0("total_", g1) := total1,
        !!paste0("n_",     g2) := n2,
        !!paste0("total_", g2) := total2,
        enrichment,
        odds_ratio,
        p_value
      )
  })
}

## ════════════════════════════════════════════════════════════════════════════
## 4.  Feature specifications & colour palette
## ════════════════════════════════════════════════════════════════════════════
features <- list(
  gene       = list(col = "gene",       lab = "Gene",        colour = c(UC="#bb9df5", CD="#ffc169")),
  product    = list(col = "product",    lab = "Product",     colour = c(UC="#bb9df5", CD="#ffc169")),
  kegg_terms = list(col = "kegg_terms", lab = "KEGG term",   colour = c(UC="#bb9df5", CD="#ffc169")),
  go_terms   = list(col = "go_terms",   lab = "GO term",     colour = c(UC="#bb9df5", CD="#ffc169")),
  ec_terms   = list(col = "ec_terms",   lab = "EC number",   colour = c(UC="#bb9df5", CD="#ffc169"))
)

## ════════════════════════════════════════════════════════════════════════════
## 6.  Run enrichment & plot volcanoes
## ════════════════════════════════════════════════════════════════════════════
## ── 0.  Put your four data-frames in a *named* list ───────────────────────
snv_sets <- list(
  snv_sweeps              = snv_sweeps,            # 1
  snv_sweeps_nosilent     = snv_sweeps_nosilent,   # 2
  snv_sweeps_1k           = snv_sweeps_1k,         # 3
  snv_sweeps_1k_nosilent  = snv_sweeps_1k_nosilent # 4
)

analyze_volcano <- function(snv_sets, features) {
  # This will become a flat named list of data.frames
  all_slices <- imap(snv_sets, function(snv_subset, sub_name) {
    # for each feature in that subset
    imap(features, function(ft, ft_name) {
      # run your enrichment
      enr <- fisher_enrich(snv_subset, ft$col)
      if (nrow(enr) == 0) return(NULL)
      
      enr2 <- enr %>%
        group_by(comparison) %>%
        mutate(
          FDR  = p.adjust(p_value, "BH"),
          logE = log10(enrichment),
          logP = -log10(p_value),
          sig  = enrichment > 2 & p_value < 0.05
        ) %>%
        ungroup()
      
      # split by comparison
      slices <- split(enr2, enr2$comparison)
      # rename each element to include subset + feature + comparison
      names(slices) <- paste(sub_name, ft_name, names(slices), sep = "_")
      slices
    }) %>%
      flatten()   # collapse the feature‐level lists
  }) %>%
    flatten()     # collapse the subset‐level lists
  
  list(data = all_slices)
}
res <- analyze_volcano(snv_sets, features)


# Now you can graph each individually
df <- res[["data"]][["snv_sweeps_1k_nosilent_gene_CD_vs_nonIBD"]]

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
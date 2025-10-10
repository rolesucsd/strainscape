## 5   Nucleotide diversity
```{r nucleotide-diversity, fig.width=2.5, fig.height=4}
select_columns <- c("Run",
                    "week_num",
                    "Participant ID",
                    "diagnosis",
                    "sex",
                    "fecalcal_ng_ml",
                    "hbi",
                    "sccai")
meta <- meta[,select_columns]

# Edit scaffold and combine with gtdbtk 
scaffold <- read.delim("/Users/reneeoles/Desktop/strainscape_output/iHMP/input/combined_processed_scaffolds.txt")
taxonomy <- read.delim("/Users/reneeoles/Desktop/strainscape_output/iHMP/input/gtdbtk.bac120.summary.tsv")
colnames(taxonomy)[1] <- "bin"
colnames(taxonomy)[21] <- "Participant ID"
colnames(scaffold)[6] <- "Run"

#Parse GTDB taxonomy into ranks
parse_gtdb <- function(df, tax_col = "classification") {
  # If your taxonomy second column has a different name, adjust above or before calling
  ranks <- c(d="d__", p="p__", c="c__", o="o__", f="f__", g="g__", s="s__")
  out <- df %>%
    mutate(
      d = str_extract(!!sym(tax_col), "d__[^;]*") %>% str_remove("^d__"),
      p = str_extract(!!sym(tax_col), "p__[^;]*") %>% str_remove("^p__"),
      c = str_extract(!!sym(tax_col), "c__[^;]*") %>% str_remove("^c__"),
      o = str_extract(!!sym(tax_col), "o__[^;]*") %>% str_remove("^o__"),
      f = str_extract(!!sym(tax_col), "f__[^;]*") %>% str_remove("^f__"),
      g = str_extract(!!sym(tax_col), "g__[^;]*") %>% str_remove("^g__"),
      s = str_extract(!!sym(tax_col), "s__[^;]*") %>% str_remove("^s__")
    )
  out
}

summary_per_sample_bin <- scaffold %>%
  group_by(Run, bin,
           Completeness, Contamination, Genome_Size) %>%
  summarise(
    n_scaffolds          = n(),
    total_length_bp      = sum(length, na.rm = TRUE),
    mean_scaffold_len_bp = mean(length, na.rm = TRUE),
    # unweighted (simple) means for reference
    coverage_mean        = mean(coverage, na.rm = TRUE),
    breadth_mean         = mean(breadth,  na.rm = TRUE),
    nucl_div_mean        = mean(nucl_diversity, na.rm = TRUE),
    coverage_max         = max(coverage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(meta) %>%
  left_join(
    taxonomy[,c(1,2,21)] %>%
      parse_gtdb(tax_col = "classification") %>%
      select(bin, classification, `Participant ID`, d, p, c, o, f, g, s)
  )

summary_per_patient <- scaffold %>%
  left_join(meta) %>%
  group_by(`Participant ID`, bin, diagnosis,
           Completeness, Contamination, Genome_Size) %>%
  summarise(
    n_scaffolds          = n(),
    total_length_bp      = sum(length, na.rm = TRUE),
    mean_scaffold_len_bp = mean(length, na.rm = TRUE),
    # unweighted (simple) means for reference
    coverage_mean        = mean(coverage, na.rm = TRUE),
    breadth_mean         = mean(breadth,  na.rm = TRUE),
    nucl_div_mean        = mean(nucl_diversity, na.rm = TRUE),
    coverage_max         = max(coverage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    taxonomy[,c(1,2,21)] %>%
      parse_gtdb(tax_col = "classification") %>%
      select(bin, classification, `Participant ID`, d, p, c, o, f, g, s)
  )
colnames(summary_per_patient)[3] <- "group"

bins_summary <- summary_per_patient %>%
  mutate(group = factor(group, levels = c("nonIBD", "UC", "CD")))

bins_summary2 <- bins_summary %>%
  mutate(
    group2 = if_else(group == "nonIBD", "nonIBD", "IBD"),
    group2 = factor(group2, levels = c("nonIBD", "IBD"))
  )

ggplot(bins_summary, aes(group, log10(nucl_div_mean), fill = group)) +
  geom_boxplot(trim = FALSE) +
  geom_hline(yintercept = c(-3, -2), linetype = "dashed", colour = "grey") +
  labs(y = "Log10 mean nucl div",
       title = "") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c("UC"     = "#bb9df5",
                               "nonIBD" = "#71d987",
                               "CD"     = "#ffc169"))
save_png(file.path(FIG_DIR, "summary/bin_nuc_div.png"), 2.5, 4)

# --- Plot log10, test on raw ---
ggplot(bins_summary2, aes(group2, log10(nucl_div_mean), fill = group2)) +
  geom_boxplot(trim = FALSE) +
  geom_hline(yintercept = c(-3, -2), linetype = "dashed", colour = "grey") +
  labs(y = "Log10 mean nucl div", title = "") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = c("IBD" = "#c474c2", "nonIBD" = "#71d987"))
save_png(file.path(FIG_DIR, "summary/bin_nuc_div_total.png"), 2, 4)

# Welch t (raw), Welch t (log10), Wilcoxon (raw),
# and attaches group-wise means/SDs (raw + log10).
pairwise_stats <- function(df, group_col, value_col, set_label = NULL, eps = NULL) {
  g <- enquo(group_col)
  y <- enquo(value_col)
  
  # epsilon only if needed (avoid log10(0))
  if (is.null(eps)) {
    miny <- min(pull(df, !!y), na.rm = TRUE)
    eps  <- if (is.finite(miny) && miny > 0) 0 else 1e-8
  }
  
  # levels/order for comparisons
  levs <- levels(pull(df, !!g))
  comps <- combn(levs, 2, simplify = FALSE)
  
  # per-group summaries
  summ <- df %>%
    group_by(!!g) %>%
    summarise(
      mean_raw = mean(!!y, na.rm = TRUE),
      sd_raw   = sd(!!y,   na.rm = TRUE),
      mean_log = mean(log10(!!y + eps), na.rm = TRUE),
      sd_log   = sd(log10(!!y + eps),   na.rm = TRUE),
      .groups = "drop"
    )
  
  # tests
  tests <- map_dfr(comps, function(comp) {
    sub <- df %>% filter(!!g %in% comp) %>% transmute(g = !!g, y = !!y)
    
    t_raw <- t.test(y ~ g,                data = sub, var.equal = FALSE) %>% tidy()
    t_log <- t.test(log10(y + eps) ~ g,   data = sub, var.equal = FALSE) %>% tidy()
    w_raw <- wilcox.test(y ~ g,           data = sub, exact = FALSE)     %>% tidy()
    
    tibble(
      group1  = comp[1],
      group2  = comp[2],
      p_t_raw = t_raw$p.value,
      p_t_log = t_log$p.value,
      p_wilcox = w_raw$p.value
    )
  }) %>%
    # attach summaries for each side
    left_join(
      summ %>% transmute(group1 = !!g,
                         mean_raw1 = mean_raw, sd_raw1 = sd_raw,
                         mean_log1 = mean_log, sd_log1 = sd_log),
      by = "group1"
    ) %>%
    left_join(
      summ %>% transmute(group2 = !!g,
                         mean_raw2 = mean_raw, sd_raw2 = sd_raw,
                         mean_log2 = mean_log, sd_log2 = sd_log),
      by = "group2"
    ) %>%
    mutate(set = set_label %||% as_label(g))
  
  tests
}

# ---- RUN BOTH SETS IN ONE GO ----

res_all <- bind_rows(
  pairwise_stats(bins_summary,  group,  nucl_div_mean, set_label = "3-group"),
  pairwise_stats(bins_summary2, group2, nucl_div_mean, set_label = "2-group (IBD vs nonIBD)")
) %>%
  group_by(set) %>%                      # adjust p-values within each family
  mutate(
    p_t_raw_adj  = p.adjust(p_t_raw,  method = "BH"),
    p_t_log_adj  = p.adjust(p_t_log,  method = "BH"),
    p_wilcox_adj = p.adjust(p_wilcox, method = "BH")
  ) %>%
  ungroup()

knitr::kable(res_all, digits = 4)



```

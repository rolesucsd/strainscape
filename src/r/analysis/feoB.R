library(data.table)
library(tidyverse)
library(scales)
library(ggridges)

bin_indentify <- read.delim("/Users/reneeoles/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Research/Strain_Evolution/iHMP/output/gtdbtk.bac120.summary.tsv")
colnames(bin_indentify)[1] <- "bin"

# all events 
ds <- open_dataset("/Users/reneeoles/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Research/Strain_Evolution/iHMP/output/all_snvs")
feoB <- ds %>% 
  filter(gene == "cirA") %>%      # row-level predicate pushed down to each file
  collect() %>%
  mutate(is_sweep = freq_range >= 0.7)
num_cols <- names(feoB)[grepl("^\\d+$", names(feoB))]
feoB_numeric <- feoB[ , c("patient_id","chromosome","position",num_cols), drop = FALSE]
feoB <- feoB[ , !(names(feoB) %in% num_cols), drop = FALSE]
feoB <- left_join(feoB, meta)
feoB <- left_join(feoB, bins_summary_filter[,c(1,2,(9:11))])
feoB <- left_join(feoB, bin_indentify[,c(1,2,21)])

table(bins_summary$group)
#CD     nonIBD     UC 
#1327    909    772 

feoB_map <- feoB %>% 
  mutate(
    is_sweep      = as.logical(is_sweep),
    mutation_type = factor(
      mutation_type,
      levels = c("Nonsense", "Missense", "Silent"),
      ordered = TRUE)
  ) %>%
  arrange(desc(is_sweep), mutation_type) %>% 
  distinct(patient_id, bin, .keep_all = TRUE)

feoB_table <- feoB_map %>% 
  mutate(mut_class = if_else(mutation_type == "Silent", "s", "ns"))
feoB_table <- as.data.frame(table(feoB_table$group, feoB_table$mut_class, feoB_table$is_sweep))
feoB_table[feoB_table$Var1 %in% "nonIBD", "Perc"] <- feoB_table[feoB_table$Var1 %in% "nonIBD", "Freq"] / 909 * 100
feoB_table[feoB_table$Var1 %in% "UC", "Perc"] <- feoB_table[feoB_table$Var1 %in% "UC", "Freq"] / 772 * 100
feoB_table[feoB_table$Var1 %in% "CD", "Perc"] <- feoB_table[feoB_table$Var1 %in% "CD", "Freq"] / 1327 * 100

feoB_table <- feoB_table %>%
  mutate(Var1 = factor(Var1, levels = c("nonIBD", "UC", "CD")))
ggplot(feoB_table, aes(x = Var2, y = Perc, fill = Var1)) +
  facet_wrap(~Var3, scales = "free_y", labeller = labeller(Var3 = c(`TRUE` = "Sweeping mutation", `FALSE` = "Non-sweeping mutation"))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("UC" = "#bb9df5",
                               "nonIBD" = "#71d987",
                               "CD" = "#ffc169")) +
  labs(x = NULL, y = "% of isolates") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 11)
  )
ggsave(file.path(OUT, "feoA_bars.png"), width = 4, height = 3.5, dpi = 320)

# test for significance
sweep_by_group <- xtabs(Freq ~ Var3 + Var1, data = feoB_table)
pairs <- combn(colnames(sweep_by_group), 2, simplify = FALSE)
pairwise_results <- lapply(pairs, function(cols) {
  tab2x2 <- sweep_by_group[ , cols]             # 2 × 2 slice
  test   <- if (all(tab2x2 >= 5))               # χ² expects ≥5 per cell
    chisq.test(tab2x2)
  else
    fisher.test(tab2x2)
  data.frame(group1 = cols[1],
             group2 = cols[2],
             p.value = test$p.value,
             test    = class(test)[1])
}) %>% bind_rows()
pairwise_results

# Ridge plot
feoB_plot <- feoB_map %>% 
  mutate(
    mut_class = if_else(mutation_type == "Silent", "s", "ns"),
    group     = factor(group, levels = c("UC", "nonIBD", "CD"))
  ) %>% 
  group_by(group, mut_class) %>%           # size of each subgroup
  mutate(w = 1 / n()) %>%                  # weight so ∑w = 1 per subgroup
  ungroup()
ggplot(
  feoB_plot,
  aes(x = freq_range, y = group, fill = group)) +
  geom_density_ridges(
    position         = "identity",         # overlay, not stacked
    scale            = 3,               # small vertical spread → more overlap
    alpha            = 0.75,              # see-through fill
    quantiles        = 2,                 # median line for each curve
    rel_min_height   = 0.01               # trim spiky tails
  ) +
  scale_fill_manual(
    values = c("UC" = "#bb9df5",
               "nonIBD" = "#71d987",
               "CD" = "#ffc169"),
    name   = "") +
  facet_wrap(~ mut_class, ncol = 1) +
  labs(x = "Frequency range",y = NULL)
ggsave(file.path(OUT, "feoA_ridge.png"), width = 3.5, height = 3.5, dpi = 320)

# BY species
feoB_species <- (as.data.frame(table(feoB_map$classification, feoB_map$group)))
feoB_species <- feoB_species %>% 
  add_count(Var1, wt = Freq, name = "total") %>%
  filter(total >= 3) %>%
  mutate(species = sub(".*s__", "", Var1)) %>%
  mutate(species = fct_reorder(species, total, .desc = FALSE))
ggplot(feoB_species, aes(fill=Var2, y=species, x=Freq)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values = c("UC" = "#bb9df5",
                               "nonIBD" = "#71d987",
                               "CD" = "#ffc169"))+
  labs(x = "Freq",y = "")
ggsave(file.path(OUT, "feoA_species.png"), width = 5, height = 5, dpi = 320)

# For each mutation that occurd intragenic, find the gene length, and the position of the mutation
feoB_mutation <- feoB %>%
  filter(coding_status == "Coding") %>%
  mutate(
    gene_length = abs(as.numeric(stop) - as.numeric(start)),
    gene_position = if_else(
      strand == "+",
      as.numeric(position) - as.numeric(start),
      as.numeric(stop) - as.numeric(position)
    )
  ) %>%
  filter(gene_length < 5000) %>%
  filter(gene_length > 1000)

ggplot(feoB_mutation, aes(x=gene_length))+
  geom_histogram() +
  facet_wrap(~gene_length) +
  scale_x_continuous(breaks = pretty_breaks(n = 20))

x_min  <- min(feoB_mutation$gene_position/3)
x_max  <- max(feoB_mutation$gene_position/3)
breaks <- seq(x_min, x_max, length.out = 101)  # 100 bins
feoB_hist <- feoB_mutation %>%
  mutate(
    x   = gene_position / 3,
    bin = cut(x, breaks = breaks, include.lowest = TRUE)
  ) %>%
  group_by(gene_length, bin) %>%
  summarise(cnt = n(), .groups = "drop")  # exactly the same as count()
keep_lengths <- feoB_hist %>%
  group_by(gene_length) %>%
  summarise(max_cnt = max(cnt), .groups="drop") %>%
  filter(max_cnt >= 5) %>%
  pull(gene_length)
feoB_filt <- feoB_mutation %>%
  filter(gene_length %in% keep_lengths)
ggplot(feoB_filt, aes(x = gene_position/3, fill = group, alpha=is_sweep)) +
  geom_histogram(breaks = breaks, color = "white") +
  #geom_vline(xintercept = c(1132, 1325, 1766, 1989, 2092)/3,linetype = "dashed") +
  scale_fill_manual(values = c(
    "UC"     = "#bb9df5",
    "nonIBD" = "#71d987",
    "CD"     = "#ffc169"
  )) +
  scale_alpha_manual(values=c(0.5,1)) +
  labs(x = NULL, y = "No. of variants") +
  facet_wrap(~ gene_length, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = pretty_breaks(8)) +
  scale_y_continuous(breaks = pretty_breaks(2))
ggsave(file.path(OUT, "udk_length_all.png"), width = 8, height = 8, dpi = 320)

feoB_mutation <- left_join(feoB_mutation, snvs_per_bin[,c(1:4)])

# Looking at the trend over time
feoB_map <- left_join(feoB, sweeps_numeric)
week_cols <- names(feoB_map)[suppressWarnings(!is.na(as.numeric(names(feoB_map))))]
feoB_map <- feoB_map %>%
  pivot_longer(
    cols = all_of(week_cols),
    names_to = "w",
    values_to = "proportion"
  )
feoB_map <- feoB_map[!is.na(feoB_map$proportion),]
feoB_map$w <- as.numeric(feoB_map$w)

feoB_map <- feoB_map %>%
  mutate(group = factor(group, levels = c("nonIBD", "UC", "CD")))
ggplot(feoB_map[feoB_map$is_sweep,], aes(
  x = w,
  y = proportion,
  group = interaction(patient_id, chromosome, position),
  color = group,
  shape = mutation_type
)) +
  geom_line()+
  geom_point()+
  scale_color_manual(values = c(
    "UC"     = "#bb9df5",
    "nonIBD" = "#71d987",
    "CD"     = "#ffc169"
  )) +
  facet_wrap(~patient_id+bin)
ggsave(file.path(OUT, "alr_sweeping.png"), width = 15, height = 15, dpi = 320)



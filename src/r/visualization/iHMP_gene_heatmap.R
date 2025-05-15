# Code for the heatmap 
combined_df_subset_metadata_summary_table <- as.data.frame(table(combined_df_subset_metadata_summary$diagnosis, combined_df_subset_metadata_summary$coding, combined_df_subset_metadata_summary$Product))
combined_df_subset_metadata_summary_table_single <- as.data.frame(table(combined_df_subset_metadata_summary$Gene))

# Convert the data into a matrix format suitable for a heatmap
heatmap_data <- combined_df_subset_metadata_summary_table[combined_df_subset_metadata_summary_table$Var2 %in% "genic",-c(2)] %>%
  select(Var1, Var3, Freq) %>%
  group_by(Var3, Var1) %>%
  summarise(Freq = sum(Freq), .groups = "drop") %>%
  pivot_wider(names_from = Var1, values_from = Freq)

# Convert to a regular data frame and set row names using Var3
heatmap_data <- as.data.frame(heatmap_data)
heatmap_data <- heatmap_data[heatmap_data$Var3 != "",]
rownames(heatmap_data) <- heatmap_data$Var3
heatmap_data <- heatmap_data %>% select(-Var3)

write.table(heatmap_data, "../iHMP/samples/wolr2/summarized_hits.txt", row.names = T, sep = "\t", quote = F)

# edited table with specific genes selected
heatmap_data <- read.delim("../iHMP/samples/wolr2/summarized_hits_edit.txt", row.names = 1)[,c(1:3)]

# Remove rows where all values are 0
heatmap_data <- heatmap_data[rowSums(heatmap_data) > 5, ]
heatmap_data_DUF <- heatmap_data[grepl("DUF", rownames(heatmap_data)), ]
heatmap_data <- heatmap_data[!grepl("tra|rib|Xer|plasmid|Trans|hypo|histidine|DNA|ecombinase|lasmid|Sigma|Phage|Mob|mob|integrase|conjugation|DUF|xerD", rownames(heatmap_data)), ]

iHMP_sample <- combined_df_subset_metadata_edit_summary_sample_species_table %>%
  count(Var2)

heatmap_data <- sweep(heatmap_data, 2, iHMP_sample$n, FUN = "/") *100

# Determine maximum value for setting breaks
max_val <- max(heatmap_data, na.rm = TRUE)
num_colors <- 50

# Create breaks so that 0 is isolated and values >=1 are scaled
breaks <- c(0, seq(1, max_val, length.out = num_colors))
# Create color palette: first color white for 0, then diverging palette for values >=1
colors <- c("white", colorRampPalette(c("#d4d3e8", "#32308a"))(num_colors - 1))

png("../iHMP/samples/wolr2/figures/summarized_product_edit.png", height = 10, width = 9, units = "in", res = 300)
# Plot the heatmap with row names, clustering, and the custom color scale
pheatmap(
  heatmap_data, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  color = colors,
  breaks = breaks,
  fontsize = 10,
  border_color = NA,
  show_rownames = TRUE
)
dev.off()


# Intergenic plot 
combined_df_subset_metadata_summary_table <- as.data.frame(table(combined_df_subset_metadata_summary$diagnosis, combined_df_subset_metadata_summary$coding, combined_df_subset_metadata_summary$Gene))
iHMP_sample <- combined_df_subset_metadata_summary_sample_species_table %>%
  count(Var2)
colnames(iHMP_sample)[1] <- "Var1"
combined_df_subset_metadata_summary_table <- left_join(combined_df_subset_metadata_summary_table, iHMP_sample)
combined_df_subset_metadata_summary_table$Freq <- combined_df_subset_metadata_summary_table$Freq/combined_df_subset_metadata_summary_table$n*100

heatmap_data <- combined_df_subset_metadata_summary_table %>%
  select(Var1, Var3, Freq) %>%
  group_by(Var3, Var1) %>%
  summarise(Freq = sum(Freq), .groups = "drop") %>%
  pivot_wider(names_from = Var1, values_from = Freq)
heatmap_data <- heatmap_data[heatmap_data$Var3 != "",]

heatmap_data <- as.data.frame(heatmap_data)
rownames(heatmap_data) <- heatmap_data$Var3
heatmap_data <- heatmap_data %>% select(-Var3)

# Now after reformatting - I am going to filter out the heatmap so that I just retain rows wehre there is sig difference between grouops
sizes <- iHMP_sample %>%
  rename(group = Var1, total_n = n)

# 2) pull heatmap_data into a tibble with a feature column
hm <- heatmap_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("feature")

# 3) reshape long, join sizes, then reshape back to wide for testing
hm2 <- hm %>%
  pivot_longer(
    cols        = c(CD, UC, nonIBD),
    names_to    = "group",
    values_to   = "pct"
  ) %>%
  left_join(sizes, by = "group") %>%    # brings in total_n for each group
  pivot_wider(
    id_cols      = feature,
    names_from   = group,
    values_from  = c(pct, total_n)
  )

# 4) run pairwise prop.tests rowwise
hm_sig <- hm2 %>%
  rowwise() %>%
  mutate(
    count_CD     = round(pct_CD     / 100 * total_n_CD),
    count_UC     = round(pct_UC     / 100 * total_n_UC),
    count_nonIBD = round(pct_nonIBD / 100 * total_n_nonIBD),
    
    p_CD_UC     = tryCatch(prop.test(
      c(count_CD, count_UC),
      c(total_n_CD, total_n_UC)
    )$p.value, error = function(e) NA_real_),
    p_CD_nonIBD = tryCatch(prop.test(
      c(count_CD, count_nonIBD),
      c(total_n_CD, total_n_nonIBD)
    )$p.value, error = function(e) NA_real_),
    p_UC_nonIBD = tryCatch(prop.test(
      c(count_UC, count_nonIBD),
      c(total_n_UC, total_n_nonIBD)
    )$p.value, error = function(e) NA_real_),
    
    # how many valid tests?
    n_valid = sum(!is.na(c(p_CD_UC, p_CD_nonIBD, p_UC_nonIBD))),
    # at least one real test and at least one p < 0.05
    any_sig = n_valid > 0 &&
      any(c(p_CD_nonIBD, p_UC_nonIBD) < 0.05,
          na.rm = TRUE)
  ) %>%
  ungroup()

stress_genes <- c(
  # master regulators / sensors
  "perR", "oxyR", "soxR", "iscR",
  
  # ROS / RNS detox
  "katE", "katG",
  "sodA", "bcp", "tpx", "rbr",
  "hmp", "dps", "btuE",
  
  # redox-recycling & Fe–S repair
  "trxB", "gor", "nifJ", "sufB", "sufC",
  
  # envelope / efflux stress
  "clpX", "clpP", "clpB", "dnaK", "dnaJ",
  "groEL", "mdtK",
  
  # iron import / export that respond to oxidative cues
  "feoB", "fieF", "fur", "fepA", "fecR", "cirA",
  
  # niche-specific candidates you mentioned earlier
  "frrB"
)

# 5) filter and pull back the original heatmap_data rows
features_keep <- hm_sig %>%
  filter(any_sig) %>%
  pull(feature)

filtered_heatmap <- heatmap_data[rownames(heatmap_data) %in% features_keep, ]
filtered_heatmap <- heatmap_data[rownames(heatmap_data) %in% stress_genes, ]

#heatmap_data <- heatmap_data[!grepl("rib|Xer|plasmid|hypo|histidine|DNA|ecombinase|lasmid|Sigma|Phage|Mob|mob|integrase|conjugation|DUF|insH|xerD|tnp|transfer|tra|trn|transposase", rownames(heatmap_data)), ]

# Determine maximum value for setting breaks
max_val <- max(filtered_heatmap, na.rm = TRUE)
num_colors <- 10

# Create breaks so that 0 is isolated and values >=1 are scaled
breaks <- c(0, seq(1, max_val, length.out = num_colors))
# Create color palette: first color white for 0, then diverging palette for values >=1
colors <- c("white", colorRampPalette(c("#d4d3e8", "#32308a"))(num_colors - 1))

png("../iHMP/samples/wolr2/figures/summarized_gene_stress_genes.png", height = 5, width = 2.5, units = "in", res = 300)
# Plot the heatmap with row names, clustering, and the custom color scale
pheatmap(
  filtered_heatmap, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  color = colors,
  breaks = breaks,
  fontsize = 10,
  border_color = NA,
  show_rownames = TRUE
)
dev.off()

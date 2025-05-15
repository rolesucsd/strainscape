combined_df_subset_metadata_genomes <- unique(combined_df_subset_metadata$Genome)

for(g in combined_df_subset_metadata_genomes) {
  df_sub <- combined_df_subset_metadata %>% filter(Genome == g)
  
  # Assumes the taxonomy string is semicolon-delimited.
  species_full <- df_sub$Species[1]  # using the first row for the genome
  species_name <- species_full %>%
    str_split(";") %>%
    unlist() %>%
    tail(1) %>%              # take the last element
    str_trim() %>%           # remove extra whitespace
    str_remove("^s__") %>%   # remove the prefix "s__" if present
    str_replace_all(" ", "_")  # replace spaces with underscores
  
  # Create the output filename.
  filename <- paste0("../iHMP/samples/wolr2/figures/species/iHMP_summary_", g, "_",species_name, ".png")
  
  # First layer: plot all points with mapping from coding and Sample.
  # Second layer: overlay points that are iron-related with a black outline.
  p <- ggplot(df_sub, aes(x = new_position, y = -log10(OLS_pvalue), starshape = coding, fill = Sample)) +
    geom_star(size = 2, alpha = 0.6, color="black") +
    scale_starshape_manual(values=c(15,13,11)) +
    facet_wrap(~diagnosis, ncol = 1) +
    custom_theme +
    #theme(legend.position = "right") +
    labs(title = "")
  
  ggsave(filename, plot = p, dpi = 300, height = 5, width = 15, units = "in")
}

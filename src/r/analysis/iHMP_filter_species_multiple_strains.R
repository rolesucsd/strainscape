#' Filter and analyze species with multiple strains
#' 
#' @description This function filters and analyzes species with multiple strains,
#' creating visualizations of strain-level data and saving filtered results.
#' 
#' @param nucleotide_div_filter Data frame containing nucleotide diversity filter
#' @param scaffold_data_edit Data frame containing scaffold data
#' @param outdir Directory to save output files
#' @param p_thresh P-value threshold for significance (default: 0.01)
#' 
#' @return Data frame containing filtered and analyzed data
#' @export
filter_species_multiple_strains <- function(nucleotide_div_filter, scaffold_data_edit, outdir, p_thresh = 0.01) {
  # Input validation
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  # Create output directories
  species_plots_dir <- file.path(outdir, "species_per_sample_plots")
  final_output_dir <- file.path(outdir, "final_output")
  dir.create(species_plots_dir, recursive = TRUE)
  dir.create(final_output_dir, recursive = TRUE)
  
  # Generate color palette
  n_colors <- 69
  random_colors <- rgb(runif(n_colors), runif(n_colors), runif(n_colors))
  df_palette <- data.frame(
    Genome = unique(nucleotide_div_filter$Genome),
    color = random_colors,
    stringsAsFactors = FALSE
  )
  named_palette <- setNames(df_palette$color, df_palette$Genome)
  
  # Process each sample
  all_large_clusters_summary <- list()
  sample_ids <- unique(nucleotide_div_filter$Sample)
  
  for(sample_id in sample_ids) {
    # Read and filter data
    df <- read_and_filter_data(sample_id, nucleotide_div_filter, scaffold_data_edit, final_output_dir)
    if (is.null(df)) next
    
    # Create visualizations
    create_species_visualizations(df, sample_id, named_palette, species_plots_dir, p_thresh)
    
    # Save filtered data
    write.table(df, 
                file.path(final_output_dir, paste0(sample_id, "_SNV_filtered_trend_trajectory_edit.txt")),
                quote = FALSE, sep = "\t", row.names = FALSE)
    
    message("Finished processing sample: ", sample_id)
  }
}

#' Read and filter data for a sample
#' @keywords internal
read_and_filter_data <- function(sample_id, nucleotide_div_filter, scaffold_data_edit, final_output_dir) {
  # Read data
  input_file <- file.path(final_output_dir, paste0(sample_id, "_SNV_filtered_trend_trajectory.txt"))
  if (!file.exists(input_file)) {
    message("Input file not found for sample: ", sample_id)
    return(NULL)
  }
  
  df <- read.delim(input_file)
  
  # Filter by nucleotide diversity and scaffold data
  allowed_species <- nucleotide_div_filter$Genome[nucleotide_div_filter$Sample == sample_id]
  allowed_chromosomes <- scaffold_data_edit$scaffold[scaffold_data_edit$Sample == sample_id]
  df <- df[df$Chromosome %in% allowed_chromosomes & df$Genome %in% allowed_species, ]
  
  if (nrow(df) == 0) {
    message("No data after filtering for sample: ", sample_id)
    return(NULL)
  }
  
  # Extract species name
  df$species <- ifelse(
    sub(".*; s__", "", df$Species) == "",
    sub(".*; g__", "", df$Species),
    sub(".*; s__", "", df$Species)
  )
  
  # Randomly shuffle the dataframe
  df[sample(nrow(df)), ]
}

#' Create visualizations for species data
#' @keywords internal
create_species_visualizations <- function(df, sample_id, named_palette, species_plots_dir, p_thresh) {
  # Calculate dimensions for faceted plot
  dim <- ceiling(sqrt(length(unique(df$species))))
  
  # Create custom labeller
  custom_labeller <- as_labeller(function(x) gsub(" ", "\n", x))
  
  # 1. Faceted plot
  p1 <- ggplot(df, aes(x = OLS_slope, y = -log10(OLS_pvalue), color=Genome)) +
    geom_point(size = 0.1) + 
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "black", alpha = 0.6, linewidth = 0.25) + 
    geom_vline(xintercept = c(-p_thresh, p_thresh), linetype = "dashed", color = "black", alpha = 0.6, linewidth = 0.25) + 
    scale_x_continuous(breaks = pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = pretty_breaks(n = 3)) +
    scale_color_manual(values = named_palette) +
    labs(x = "OLS Slope", y = expression(-log[10]("P-value"))) +
    theme_ihmp() +  
    facet_wrap(~species, labeller = custom_labeller)
  
  ggsave(file.path(species_plots_dir, paste0(sample_id, "_facet.png")),
         plot = p1,
         dpi = 300, height = dim*2.5, width = dim*2.5, units = "in")
  
  # 2. Combined plot
  p2 <- ggplot(df, aes(x = OLS_slope, y = -log10(OLS_pvalue), color = Genome)) +
    geom_jitter(aes(
      size = (-log10(OLS_pvalue) > -log10(p_thresh)),
      alpha = (-log10(OLS_pvalue) > -log10(p_thresh))
    ), width = 0, height = 0.1) +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-p_thresh, p_thresh), linetype = "dashed", color = "black") + 
    scale_color_manual(values = named_palette) +
    labs(x = "OLS Slope", y = expression(-log[10]("P-value"))) +
    scale_size_manual(values = c(0.2, 0.5)) +
    scale_alpha_manual(values = c(0.5, 1)) +
    theme_ihmp()
  
  ggsave(file.path(species_plots_dir, paste0(sample_id, "_combined.png")),
         plot = p2,
         dpi = 300, height = 7, width = 7, units = "in")
} 
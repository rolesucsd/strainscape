#!/usr/bin/env Rscript

# iHMP Analysis Pipeline - R Analysis Wrapper
# This script runs all R analyses on the processed data

suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(tidyverse)
  library(data.table)
  library(progress)
  library(ggpubr)
  library(patchwork)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--config", "-c"), type="character", default="config.yaml",
              help="Path to configuration file"),
  make_option(c("--input", "-i"), type="character", default="results/intermediate",
              help="Input directory containing processed data"),
  make_option(c("--output", "-o"), type="character", default="results/final",
              help="Output directory for analysis results")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Source analysis functions
source(here("src", "r", "analysis", "analyze_scaffolds.R"))
source(here("src", "r", "analysis", "analyze_mutations.R"))
source(here("src", "r", "analysis", "analyze_genes.R"))
source(here("src", "r", "visualization", "create_plots.R"))

# Create output directories
dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$output, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$output, "tables"), recursive = TRUE, showWarnings = FALSE)

# Run analyses
cat("Running scaffold analysis...\n")
scaffold_results <- analyze_scaffolds(
  input_dir = opt$input,
  output_dir = opt$output
)

cat("Running mutation analysis...\n")
mutation_results <- analyze_mutations(
  input_dir = opt$input,
  output_dir = opt$output
)

cat("Running gene analysis...\n")
gene_results <- analyze_genes(
  input_dir = opt$input,
  output_dir = opt$output
)

cat("Creating visualizations...\n")
create_all_plots(
  scaffold_results = scaffold_results,
  mutation_results = mutation_results,
  gene_results = gene_results,
  output_dir = file.path(opt$output, "figures")
)

cat("Analysis complete!\n")
cat("Results saved to:", opt$output, "\n") 
# -----------------------------
# iHMP Genome Summary Analysis
# -----------------------------
# This script processes and summarizes genome-wide mutation data from iHMP:
# 1. Filters significant mutations based on p-value and effect size
# 2. Joins with scaffold metadata and sample information
# 3. Expands database cross-references for functional annotation
# 4. Generates summary statistics and visualizations

# Libraries
library(tidyverse)
library(data.table)    # for fast I/O and operations
library(progress)      # for progress bars

# ----- Configuration & Paths -----
INPUT_FN   <- "combined_df.txt"
REF_SUM_FN <- file.path("..", "iHMP", "samples", "wolr2", "scaffold_info", "nucleotide_div_filter.txt")
OUTPUT_SUM <- file.path("..", "iHMP", "samples", "wolr2", "summarized_hits_by_species_and_sample.txt")
META_FP    <- file.path("..", "iHMP", "metadata", "hmp2_metadata_2018-08-20.csv")

# ----- Constants with Biological Significance -----
# SIG_P: Significance threshold for p-values (0.05)
# - Standard threshold for statistical significance
# - Balances false positives with biological relevance
SIG_P      <- 0.05

# SIG_SLOPE: Minimum effect size threshold (0.001)
# - Filters out mutations with minimal impact
# - Ensures biological relevance of detected changes
SIG_SLOPE  <- 0.001

# ----- Helper Functions -----
read_and_filter_combined_data <- function(input_fn) {
  # Read and filter combined data with progress bar
  pb <- progress_bar$new(total = 3, format = "Reading combined data [:bar] :percent")
  pb$tick()
  
  # Read and initial filtering
  combined_df <- fread(input_fn, sep = "\t", colClasses = "character") %>%
    # drop unneeded columns
    .[, !c("V10", "V11", "V12", "10", "11", "12"), with = FALSE] %>%
    # filter raw p-values
    .[, OLS_pvalue := as.numeric(OLS_pvalue)] %>%
    .[OLS_pvalue <= 0.10]
  
  pb$tick()
  
  # Join scaffold metadata
  ref_meta <- fread(REF_SUM_FN, sep = "\t", colClasses = "character") %>%
    setnames(c("V7", "V1"), c("Chromosome", "Sample")) %>%
    .[, Sample := str_remove(Sample, "\\.txt$")]
  
  combined_df <- combined_df[ref_meta, on = "Sample", nomatch = 0]
  
  pb$tick()
  
  combined_df
}

process_mutation_data <- function(combined_df, meta_fp) {
  # Process mutation data with progress bar
  pb <- progress_bar$new(total = 3, format = "Processing mutation data [:bar] :percent")
  pb$tick()
  
  # Add SNP IDs and process metadata
  combined_df <- combined_df %>%
    .[, `:=`(
      snp_id = paste(Chromosome, Position, sep = "_"),
      Mutation_Type = fcoalesce(Mutation_Type, "Intergenic"),
      diagnosis = factor(diagnosis, levels = c("nonIBD", "CD", "UC")),
      Position = as.numeric(Position)
    )]
  
  pb$tick()
  
  # Join sample metadata
  meta <- fread(meta_fp) %>%
    .[, External.ID := str_remove(External.ID, "_P$")] %>%
    .[, .(Sample = External.ID, diagnosis, week)]
  
  combined_df <- combined_df[meta, on = "Sample"]
  
  pb$tick()
  
  combined_df
}

summarize_significant_hits <- function(combined_df) {
  # Summarize significant hits with progress bar
  pb <- progress_bar$new(total = 3, format = "Summarizing hits [:bar] :percent")
  pb$tick()
  
  # Filter significant hits
  summary_df <- combined_df[
    OLS_pvalue <= SIG_P & abs(OLS_slope) > SIG_SLOPE,
    .(Sample, Chromosome, Position, upstream_Product, downstream_Product,
      Matched_Product, coding)
  ] %>%
    unique()
  
  pb$tick()
  
  # Pivot to long format
  summary_df <- summary_df %>%
    melt(
      measure.vars = patterns("_Product$", "_Gene$", "_DbXrefs$"),
      variable.name = "Location",
      value.name = c("Product", "Gene", "DbXrefs")
    ) %>%
    .[Product != ""]
  
  pb$tick()
  
  summary_df
}

expand_dbxrefs <- function(summary_df) {
  # Expand database cross-references with progress bar
  pb <- progress_bar$new(total = 2, format = "Expanding cross-references [:bar] :percent")
  pb$tick()
  
  # Split and expand DbXrefs
  expanded_df <- summary_df %>%
    .[, .(DbXrefs = unlist(strsplit(DbXrefs, ",\\s*"))), by = setdiff(names(summary_df), "DbXrefs")] %>%
    .[DbXrefs != ""] %>%
    .[, c("key", "value") := tstrsplit(DbXrefs, ":", fixed = TRUE)] %>%
    .[!is.na(key) & key != "" & !is.na(value) & value != ""] %>%
    dcast(
      ... ~ key,
      value.var = "value",
      fun.aggregate = function(x) paste(unique(x), collapse = ";"),
      fill = NA_character_
    )
  
  pb$tick()
  
  expanded_df
}

# ----- Main Analysis -----
# Read and process data
combined_data <- read_and_filter_combined_data(INPUT_FN)
mutation_data <- process_mutation_data(combined_data, META_FP)
significant_hits <- summarize_significant_hits(mutation_data)
expanded_hits <- expand_dbxrefs(significant_hits)

# ----- Write Output -----
fwrite(expanded_hits, OUTPUT_SUM, sep = "\t")

# ----- Generate Summary Statistics -----
# Count table of mutations by sample and species
count_tbl <- expanded_hits[
  , .N,
  by = .(Sample, diagnosis, Species, Chromosome, Position)
][N > 0]

# Optional: filter high-count combinations
# MAX_MUT <- 1000
# overcount <- count_tbl[N >= MAX_MUT]

# ----- Visualization (if needed) -----
# Add visualization code here if required
# Example:
# count_tbl %>%
#   ggplot(aes(diagnosis, N, fill = diagnosis)) +
#   geom_boxplot() +
#   facet_wrap(~ Species) +
#   labs(title = "Mutation Count Distribution by Diagnosis",
#        subtitle = "Number of significant mutations per sample",
#        x = "Diagnosis",
#        y = "Number of Mutations") +
#   custom_theme

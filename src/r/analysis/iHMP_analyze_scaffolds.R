# -----------------------------
# iHMP Scaffold & Diversity Analysis
# -----------------------------
# This script analyzes metagenomic data from the iHMP dataset, focusing on:
# 1. Scaffold quality and coverage analysis
# 2. Nucleotide diversity patterns across different diagnoses (UC, CD, nonIBD)
# 3. Temporal dynamics of microbial diversity

# Libraries
library(tidyverse)
library(data.table)    # for fast I/O and operations
library(progress)      # for progress bars
library(ggpubr)        # for statistical tests
library(patchwork)     # for plot arrangement

# ----- Configuration & Paths -----
ROOT    <- "../iHMP"
META_FP <- file.path(ROOT, "metadata", "hmp2_metadata_2018-08-20.csv")
SPECIES_LIST_FP <- file.path(ROOT, "samples", "wolr2", "reference", "species_list.txt")
FLAGSTAT_FP     <- file.path(ROOT, "samples", "wolr2", "scaffold_info", "combined.flagstat.txt")
SCAFFOLD_FP     <- file.path(ROOT, "samples", "wolr2", "scaffold_info", "combined_scaffold_info.tsv")
FIG_DIR <- file.path(ROOT, "samples", "wolr2", "figures", "scaffolds")

# ----- Constants with Biological Significance -----
# MIN_SCAF_LEN: Minimum scaffold length (1000bp)
# - Shorter scaffolds may represent assembly artifacts or low-quality regions
# - 1000bp threshold ensures sufficient sequence for reliable diversity estimation
MIN_SCAF_LEN    <- 1000      

# MIN_WEEK_LEN: Minimum coverage per week (1Mb)
# - Ensures sufficient sequencing depth for reliable diversity estimates
# - 1Mb coverage provides good statistical power for nucleotide diversity calculations
MIN_WEEK_LEN    <- 1e6       

# MIN_GOOD_WEEKS: Minimum number of weeks with good coverage (3)
# - Ensures temporal stability of the genome
# - 3 weeks provides sufficient data points for temporal trend analysis
MIN_GOOD_WEEKS  <- 3         

# ----- Helper Functions -----
read_and_clean_metadata <- function(meta_fp) {
  meta <- fread(meta_fp) %>%
    # standardize Sample IDs
    .[, External.ID := str_remove(External.ID, "_P$")] %>%
    .[data_type == "metagenomics"]  # keep only metagenomics samples
  
  # count samples per participant
  sample_counts <- meta[, .N, by = `Participant ID`][N >= 5]
  
  # pull simple lookup for diagnosis
  meta_diag <- unique(meta[, .(`Participant ID`, diagnosis)])
  
  list(metadata = meta, sample_counts = sample_counts, meta_diag = meta_diag)
}

process_flagstat <- function(flagstat_fp, sample_counts) {
  flagstat <- fread(flagstat_fp, header = FALSE,
                   col.names = c("MappedPct","ID","Sample")) %>%
    .[, Mapped := parse_number(MappedPct)] %>%
    .[Sample %in% sample_counts$`Participant ID`]
  
  flagstat
}

process_scaffold_data <- function(scaffold_fp, species_list_fp, metadata, sample_counts) {
  # Read and clean scaffold data
  scaffold <- fread(scaffold_fp, sep = "\t", check.names = FALSE, fill = TRUE) %>%
    # Clean column names
    setnames(names(.), gsub("\\s+", "_", names(.))) %>%
    # Remove any empty rows
    .[!is.na(scaffold)] %>%
    # Clean up scaffold names
    .[, scaffold := gsub("\\s+", "", scaffold)] %>%
    # Convert length to numeric
    .[, length := as.numeric(gsub("\\s+", "", length))] %>%
    # Filter by minimum length
    .[length > MIN_SCAF_LEN] %>%
    # Extract genome ID
    .[, Genome := str_replace(Assembly, "^GC[FA]_(\\d+)\\.[12]$", "G\\1")]
  
  scaffold <- scaffold[,-c(25:27)]
  
  # Attach taxonomy with correct column names
  species_list <- fread(species_list_fp, 
                       col.names = c("Genome", "K", "P", "C", "O", "F", "G", "S", "S2"))
  scaffold <- scaffold[species_list, on = "Genome"]
  
  # Clean up sample IDs and attach metadata
  scaffold <- scaffold %>%
    .[, Sample := gsub("\\s+", "", Sample)] %>%
    .[metadata, on = c("Sample" = "External.ID")] %>%
    .[`Participant ID` %in% sample_counts$`Participant ID`]
  
  scaffold
}

filter_by_coverage <- function(scaffold) {
  # Compute total scaffold length per (Participant, Genome, week)
  week_cov <- scaffold[, .(week_length = sum(length, na.rm = TRUE)), 
                      by = .(`Participant ID`, Genome, week_num)]
  
  # Find good genomes
  good_genomes <- week_cov[week_length >= MIN_WEEK_LEN] %>%
    .[, .(n_weeks = uniqueN(week_num)), by = .(`Participant ID`, Genome)] %>%
    .[n_weeks >= MIN_GOOD_WEEKS] %>%
    .[, .(`Participant ID`, Genome)]

    # Filter scaffold data
  scaffold_filt <- scaffold[good_genomes, on = c("Participant ID", "Genome")] %>%
    .[, Species := str_extract(paste(S, S2), "(?<=s__).*")] %>%
    .[Species == "", Species := NA]
  
  scaffold_filt
}

analyze_diversity <- function(scaffold_filt, species = "ovatus") {
  # Per-species summary with exact matching
  species_df <- scaffold_filt[S2 == species] %>%
    .[, .(mean_div = mean(nucl_diversity, na.rm = TRUE),
          se_div = sd(nucl_diversity, na.rm = TRUE)/sqrt(.N),
          n_pts = .N),
      by = .(diagnosis, `Participant ID`, S2, week_num)]
  
  # First filter to keep only participants with enough total timepoints
  valid_participants <- species_df[, .(total_pts = .N), by = .(`Participant ID`)] %>%
    .[total_pts >= 3, `Participant ID`]
  
  # Then filter the main dataset to keep only valid participants
  species_df <- species_df[`Participant ID` %in% valid_participants]
  
  # Low diversity analysis
  low_div <- scaffold_filt[, .(med_div = median(nucl_diversity, na.rm = TRUE),
                              se_div = sd(nucl_diversity, na.rm = TRUE)/sqrt(.N)),
                          by = .(`Participant ID`, Genome, S2, diagnosis)] %>%
    .[med_div < 0.001 & !is.na(Genome)]
  
  list(species_df = species_df, low_div = low_div)
}

# ----- Main Analysis -----
# Read and process data
meta_data <- read_and_clean_metadata(META_FP)
flagstat_data <- process_flagstat(FLAGSTAT_FP, meta_data$sample_counts)
scaffold_data <- process_scaffold_data(SCAFFOLD_FP, SPECIES_LIST_FP, 
                                     meta_data$metadata, meta_data$sample_counts)
scaffold_filt <- filter_by_coverage(scaffold_data)
diversity_data <- analyze_diversity(scaffold_filt)

# ----- Enhanced Plots -----
# Flagstat QC Plot with enhanced titles
flagstat_data %>%
  ggplot(aes(Mapped)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ Sample) +
  labs(title = "Mapping Quality Distribution by Sample",
       subtitle = "Distribution of percentage of reads mapped to reference genomes",
       x = "Percentage of Reads Mapped",
       y = "Count") +
  custom_theme

# Nucleotide diversity over time with enhanced titles and statistical tests
diversity_data$species_df %>%
  mutate(`Participant ID` = factor(`Participant ID`, 
                                 levels = unique(`Participant ID`))) %>%
  ggplot(aes(week_num, mean_div*100, color = diagnosis, group = `Participant ID`)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = (mean_div - se_div)*100,
                    ymax = (mean_div + se_div)*100),
                width = 0.3) +
  facet_wrap(~ `Participant ID`, scales = "free_y") +
  scale_color_manual(values = c("UC"="#bb9df5",
                               "nonIBD"="#71d987",
                               "CD"="#ffc169")) +
  labs(title = "Temporal Dynamics of Nucleotide Diversity",
       subtitle = "Bacteroides ovatus diversity patterns across different diagnoses",
       x = "Week",
       y = "Nucleotide Diversity (%)") +
  custom_theme

# Enhanced statistical tests for diversity comparisons
div_median <- scaffold_filt[, .(mean_div = median(nucl_diversity, na.rm = TRUE)),
                           by = .(diagnosis, `Participant ID`, Species, Genome)] %>%
  .[!is.na(Species)] %>%
  .[, Species := str_replace_all(Species, " ", "\n")]

# Save plots
ggsave(file.path(FIG_DIR, "nucl_diversity_ovatus.png"),
       height = 12, width = 20, units = "in")

# Perform multiple statistical tests
valid_species <- div_median[, .N, by = .(Species, diagnosis)] %>%
  .[N >= 3, unique(Species)]

# Create filtered plot
div_median[Species %in% valid_species] %>%
  ggplot(aes(diagnosis, log10(mean_div), fill = diagnosis)) +
  geom_violin() + 
  geom_jitter(width = 0.1, alpha=0.5) +
  geom_hline(yintercept = c(-3,-2), linetype = "dashed", color = "grey") +
  stat_compare_means(comparisons = list(c("nonIBD","UC"),
                                       c("nonIBD","CD"),
                                       c("UC","CD")),
                    method = "wilcox.test", label = "p.format") +
  stat_compare_means(method = "kruskal.test", label.y = -0.5) +  # Add Kruskal-Wallis test
  facet_wrap(~ Species, scales = "free_y") +
  coord_cartesian(ylim = c(-3.5,-0.5)) +
  scale_fill_manual(values = c("UC"="#bb9df5",
                              "nonIBD"="#71d987",
                              "CD"="#ffc169")) +
  labs(title = "Nucleotide Diversity Distribution by Diagnosis",
       subtitle = "Violin plots showing distribution of median nucleotide diversity across diagnoses",
       x = "Diagnosis",
       y = "Log10(Nucleotide Diversity)") +
  custom_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(FIG_DIR, "nucl_diversity_median_summarized_species.png"),
       height=16, width=12.5, units="in")


div_median %>%
  ggplot(aes(diagnosis, log10(mean_div), fill = diagnosis)) +
  geom_boxplot() + 
  geom_hline(yintercept = c(-3,-2), linetype = "dashed", color = "grey") +
  stat_compare_means(comparisons = list(c("nonIBD","UC"),
                                        c("nonIBD","CD"),
                                        c("UC","CD")),
                     method = "wilcox.test", label = "p.format") +
  stat_compare_means(method = "kruskal.test", label.y = -0.5) +  # Add Kruskal-Wallis test
  coord_cartesian(ylim = c(-3.5,-0.5)) +
  scale_fill_manual(values = c("UC"="#bb9df5",
                               "nonIBD"="#71d987",
                               "CD"="#ffc169")) +
  labs(title = "Nucleotide Diversity Distribution by Diagnosis",
       subtitle = "Violin plots showing distribution of median nucleotide diversity across diagnoses",
       x = "Diagnosis",
       y = "Log10(Nucleotide Diversity)") +
  custom_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(FIG_DIR, "nucl_diversity.png"),
       height=4, width=3, units="in")

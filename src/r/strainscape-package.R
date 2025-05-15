#' StrainScape: Strain Evolution Analysis and Visualization
#'
#' A package for analyzing and visualizing strain evolution in microbiome data.
#' Provides tools for processing, analyzing, and visualizing strain-level changes
#' in microbial communities.
#'
#' @docType package
#' @name strainscape
#' @importFrom magrittr %>%
#' @importFrom dplyr select filter mutate group_by summarize
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme_minimal
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom readr read_csv write_csv
#' @importFrom stringr str_detect str_replace
NULL 
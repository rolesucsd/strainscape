#' Load Data
#' 
#' @description Loads data from various file formats.
#' 
#' @param file_path Path to the data file
#' @param file_type Type of file (default: "csv")
#' 
#' @return A data frame containing the loaded data
#' 
#' @export
load_data <- function(file_path, file_type = "csv") {
  if (!file.exists(file_path)) {
    stop(sprintf("File not found: %s", file_path))
  }
  
  switch(file_type,
         "csv" = read.csv(file_path, stringsAsFactors = FALSE),
         "tsv" = read.delim(file_path, stringsAsFactors = FALSE),
         "rds" = readRDS(file_path),
         stop(sprintf("Unsupported file type: %s", file_type)))
}

#' Save Data
#' 
#' @description Saves data to various file formats.
#' 
#' @param data Data to save
#' @param file_path Path to save the data
#' @param file_type Type of file (default: "csv")
#' 
#' @export
save_data <- function(data, file_path, file_type = "csv") {
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  
  switch(file_type,
         "csv" = write.csv(data, file_path, row.names = FALSE),
         "tsv" = write.table(data, file_path, sep = "\t", row.names = FALSE),
         "rds" = saveRDS(data, file_path),
         stop(sprintf("Unsupported file type: %s", file_type)))
}

#' Filter Data
#' 
#' @description Filters data based on specified criteria.
#' 
#' @param data Data frame to filter
#' @param filters List of filtering criteria
#' 
#' @return Filtered data frame
#' 
#' @export
filter_data <- function(data, filters) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  if (!is.list(filters)) {
    stop("filters must be a list")
  }
  
  filtered_data <- data
  for (filter_name in names(filters)) {
    if (!filter_name %in% colnames(filtered_data)) {
      warning(sprintf("Column not found: %s", filter_name))
      next
    }
    
    filter_value <- filters[[filter_name]]
    if (is.numeric(filter_value)) {
      filtered_data <- filtered_data[filtered_data[[filter_name]] >= filter_value, ]
    } else {
      filtered_data <- filtered_data[filtered_data[[filter_name]] %in% filter_value, ]
    }
  }
  
  filtered_data
}

#' Merge Data
#' 
#' @description Merges multiple data frames.
#' 
#' @param data_list List of data frames to merge
#' @param by Column(s) to merge by
#' @param all Whether to keep all rows (default: TRUE)
#' 
#' @return Merged data frame
#' 
#' @export
merge_data <- function(data_list, by, all = TRUE) {
  if (!is.list(data_list)) {
    stop("data_list must be a list")
  }
  if (length(data_list) < 2) {
    stop("data_list must contain at least 2 data frames")
  }
  
  merged_data <- data_list[[1]]
  for (i in 2:length(data_list)) {
    merged_data <- merge(merged_data, data_list[[i]], by = by, all = all)
  }
  
  merged_data
} 
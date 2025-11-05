# =====================================
# Task 1 Functions
# =====================================

#' Load count and metadata files
#'
#' Reads two input files: a counts table and a metadata table.
#'
#' @param counts_file Path to the counts file (CSV or TSV) containing expression counts.
#' @param metadata_file Path to the metadata file (CSV or TSV) containing sample information.
#' @return A list with two elements: counts and metadata as data.tables.
#' @export
load_data <- function(counts_file, metadata_file) {
  counts <- data.table::fread(counts_file)
  metadata <- data.table::fread(metadata_file)
  list(counts = counts, metadata = metadata)
}

#' Filter counts and metadata based on gene pattern and condition
#'
#' @param counts A data.table containing gene count data with a `gene` column.
#' @param metadata A data.table containing sample metadata with a `condition` column.
#' @param gene_pattern A regular expression pattern to match gene names.
#' @param condition_value A character string indicating the condition to filter samples by.
#' @return A list with filtered `counts` and `metadata`.
#' @export
filter_data <- function(counts, metadata, gene_pattern, condition_value) {
  counts_filtered <- counts[gene %like% gene_pattern]
  metadata_filtered <- metadata[condition == condition_value]
  list(counts = counts_filtered, metadata = metadata_filtered)
}

#' Merge filtered counts and metadata by sample_id
#'
#' @param counts_filtered Filtered counts data.table with `sample_id`.
#' @param metadata_filtered Filtered metadata data.table with `sample_id`.
#' @return A merged data.table.
#' @export
merge_data <- function(counts_filtered, metadata_filtered) {
  data.table::merge(counts_filtered, metadata_filtered, by = "sample_id")
}

#' Summarize merged data by gene
#'
#' @param DT A merged data.table containing `gene` and `count` columns.
#' @return A data.table with mean and median counts per gene.
#' @export
summarize_data <- function(DT) {
  DT[, .(
    mean_count = mean(count),
    median_count = median(count)
  ), by = gene]
}

#' Print a summary of top genes
#'
#' @param summary_dt A data.table returned by `summarize_data()`.
#' @param n Number of rows to display. Defaults to 5.
#' @return Invisibly returns the printed subset.
#' @export
print_summary <- function(summary_dt, n = 5) {
  print(head(summary_dt, n))
  invisible(head(summary_dt, n))
}

#' Run Task 1: Load, filter, merge, summarize, and print data
#'
#' @param counts_file Path to the counts file. Defaults to `"bulk_counts_long.csv"`.
#' @param metadata_file Path to the metadata file. Defaults to `"sample_metadata.csv"`.
#' @param gene_pattern Regex to match gene names. Defaults to `"^GENE_00"`.
#' @param condition_value Condition to filter metadata by. Defaults to `"treated"`.
#' @return Invisibly returns the summary data.table.
#' @export
task1 <- function(counts_file = "bulk_counts_long.csv",
                  metadata_file = "sample_metadata.csv",
                  gene_pattern = "^GENE_00",
                  condition_value = "treated") {
  data <- load_data(counts_file, metadata_file)
  filtered <- filter_data(data$counts, data$metadata, gene_pattern, condition_value)
  merged <- merge_data(filtered$counts, filtered$metadata)
  summary_dt <- summarize_data(merged)
  print_summary(summary_dt)
  invisible(summary_dt)
}

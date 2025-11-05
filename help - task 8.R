#' Load bulk RNA-seq counts and metadata
#'
#' Reads count and metadata tables from CSV files and returns them as a list of `data.table`s.
#'
#' @param counts_file Character string. Path to the counts file (default = `"bulk_counts_long.csv"`).
#'   Expected columns: `gene`, `sample_id`, `count`.
#' @param metadata_file Character string. Path to the metadata file (default = `"sample_metadata.csv"`).
#'   Expected columns: `sample_id`, `condition`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{counts}{`data.table` containing count data.}
#'   \item{metadata}{`data.table` containing sample metadata.}
#' }
#'
#' @export
load_task8_data <- function(counts_file = "bulk_counts_long.csv",
                            metadata_file = "sample_metadata.csv") {
  counts <- fread(counts_file)
  metadata <- fread(metadata_file)
  return(list(counts = counts, metadata = metadata))
}

#' Merge count and metadata tables
#'
#' Merges bulk RNA-seq count data with sample metadata by the `sample_id` column.
#'
#' @param counts A `data.table` containing count data.
#' @param metadata A `data.table` containing sample metadata.
#'
#' @return A merged `data.table` with combined count and metadata information.
#'
#' @export
merge_counts_metadata <- function(counts, metadata) {
  dt <- merge(counts, metadata, by = "sample_id")
  return(dt)
}

#' Compute per-gene summary statistics by condition
#'
#' Calculates summary statistics for each gene across conditions.
#'
#' @param dt A merged `data.table` containing columns: `gene`, `condition`, and `count`.
#'
#' @return A `data.table` containing, for each gene and condition:
#' \describe{
#'   \item{mean_count}{Mean count.}
#'   \item{median_count}{Median count.}
#'   \item{Q1}{First quartile (25th percentile).}
#'   \item{Q3}{Third quartile (75th percentile).}
#' }
#'
#' @export
compute_summary_stats <- function(dt) {
  summary_dt <- dt[, .(
    mean_count = mean(count, na.rm = TRUE),
    median_count = median(count, na.rm = TRUE),
    Q1 = quantile(count, 0.25, na.rm = TRUE),
    Q3 = quantile(count, 0.75, na.rm = TRUE)
  ), by = .(gene, condition)]
  return(summary_dt)
}

#' Reshape summary statistics into wide format
#'
#' Converts a long-format summary table into a wide format, with one column per condition.
#'
#' @param summary_dt A `data.table` produced by `compute_summary_stats()`.
#' @param value_col Character string. Column to use as the value in the wide table (default = `"mean_count"`).
#'
#' @return A wide-format `data.table` with one row per gene and one column per condition.
#'
#' @export
reshape_summary_wide <- function(summary_dt, value_col = "mean_count") {
  summary_wide <- dcast(summary_dt, gene ~ condition, value.var = value_col)
  return(summary_wide)
}

#' Filter genes based on fold change between conditions
#'
#' Identifies genes whose expression in the treated condition is at least
#' a specified multiple of the control condition.
#'
#' @param summary_wide A wide-format `data.table` containing mean counts by condition.
#' @param control_name Character string. Name of the control condition column (default = `"control"`).
#' @param treated_name Character string. Name of the treated condition column (default = `"treated"`).
#' @param fold_change Numeric. Minimum fold change threshold (default = 2).
#'
#' @return A filtered `data.table` containing genes meeting the fold-change criterion,
#'   along with control and treated mean counts.
#'
#' @export
filter_genes_fc <- function(summary_wide, control_name = "control", treated_name = "treated", fold_change = 2) {
  filtered_genes <- summary_wide[get(treated_name) >= fold_change * get(control_name),
                                 .(gene, control = get(control_name), treated = get(treated_name))]
  return(filtered_genes)
}

#' Merge filtered genes with additional summary statistics
#'
#' Adds median and quartile values for the filtered genes.
#'
#' @param filtered_genes A `data.table` produced by `filter_genes_fc()`.
#' @param summary_dt A long-format summary `data.table` produced by `compute_summary_stats()`.
#'
#' @return A merged `data.table` containing filtered genes along with their median, Q1, and Q3 per condition.
#'
#' @export
merge_other_stats <- function(filtered_genes, summary_dt) {
  other_stats <- dcast(summary_dt, gene ~ condition, value.var = c("median_count", "Q1", "Q3"))
  result <- merge(filtered_genes, other_stats, by = "gene")
  return(result)
}

#' Bulk RNA-seq differential expression summary pipeline
#'
#' Computes per-gene summary statistics, identifies genes with fold-change differences,
#' merges with additional statistics, and saves the results.
#'
#' @param counts_file Character string. Path to the counts file (default = `"bulk_counts_long.csv"`).
#' @param metadata_file Character string. Path to the metadata file (default = `"sample_metadata.csv"`).
#' @param output_file Character string. Path to save the filtered summary results (default = `"filtered_gene_summary.csv"`).
#'
#' @return Invisibly returns the result `data.table` with filtered genes and summary statistics.
#'
#' @details
#' The pipeline performs the following steps:
#' \enumerate{
#'   \item Load count and metadata tables.
#'   \item Merge them by `sample_id`.
#'   \item Compute mean, median, and quantile statistics per gene per condition.
#'   \item Reshape to wide format and apply fold-change filtering.
#'   \item Merge additional statistics and export results.
#' }
#'
#' @export
task8 <- function(counts_file = "bulk_counts_long.csv",
                  metadata_file = "sample_metadata.csv",
                  output_file = "filtered_gene_summary.csv") {
  
  data_list <- load_task8_data(counts_file, metadata_file)
  dt <- merge_counts_metadata(data_list$counts, data_list$metadata)
  
  summary_dt <- compute_summary_stats(dt)
  summary_wide <- reshape_summary_wide(summary_dt, value_col = "mean_count")
  
  filtered_genes <- filter_genes_fc(summary_wide)
  result <- merge_other_stats(filtered_genes, summary_dt)
  
  fwrite(result, output_file)
  print(head(result))
  
  invisible(result)
}

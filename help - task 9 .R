#' Load bulk RNA-seq counts in wide format
#'
#' Reads a CSV file containing bulk RNA-seq count data in wide format.
#'
#' @param filepath Character string. Path to the counts file (default = `"bulk_counts_wide.csv"`).
#'
#' @return A `data.table` containing wide-format counts, where each row corresponds to a gene
#'   and each column (except the ID) corresponds to a sample.
#'
#' @export
load_wide_counts <- function(filepath = "bulk_counts_wide.csv") {
  bulk_counts <- fread(filepath)
  return(bulk_counts)
}

#' Convert wide-format counts to long format
#'
#' Converts a wide-format count matrix to long format, suitable for aggregation and normalization.
#'
#' @param bulk_counts A `data.table` containing wide-format counts.
#' @param id_col Character string. Column name identifying the unique ID variable (default = `"gene"`).
#' @param value_name Character string. Name of the value column in the long format (default = `"count"`).
#'
#' @return A long-format `data.table` with columns: `gene`, `sample`, and `count`.
#'
#' @export
convert_wide_to_long <- function(bulk_counts, id_col = "gene", value_name = "count") {
  long_counts <- melt(
    bulk_counts,
    id.vars = id_col,
    variable.name = "sample",
    value.name = value_name
  )
  return(long_counts)
}

#' Normalize counts by total per sample
#'
#' Computes the total count per sample and adds normalized counts to the data.
#'
#' @param long_counts A long-format `data.table` containing columns `sample` and `count`.
#'
#' @return A `data.table` with added columns:
#' \describe{
#'   \item{total_count}{Total counts per sample.}
#'   \item{normalized_count}{Count normalized by the sample total.}
#' }
#'
#' @export
add_sample_totals <- function(long_counts) {
  sample_totals <- long_counts[, .(total_count = sum(count)), by = sample]
  long_counts <- merge(long_counts, sample_totals, by = "sample")
  long_counts[, normalized_count := count / total_count]
  return(long_counts)
}

#' Extract experimental condition from sample names
#'
#' Splits the sample identifier string to extract a condition label.
#'
#' @param long_counts A `data.table` containing a `sample` column.
#' @param pattern_split Character string used to split the sample names (default = `"_"`).
#' @param part Integer specifying which part of the split name to use (default = `1`).
#'
#' @return A `data.table` with an added column `condition`.
#'
#' @export
extract_condition <- function(long_counts, pattern_split = "_", part = 1) {
  long_counts[, condition := tstrsplit(sample, pattern_split)[[part]]]
  return(long_counts)
}

#' Compute mean counts per gene and condition
#'
#' Aggregates normalized (or raw) counts to calculate the mean count for each gene and condition.
#'
#' @param long_counts A `data.table` containing columns `gene`, `condition`, and `count`.
#'
#' @return A `data.table` with columns `gene`, `condition`, and `mean_count`.
#'
#' @export
compute_gene_condition_means <- function(long_counts) {
  gene_condition_means <- long_counts[, .(mean_count = mean(count, na.rm = TRUE)),
                                      by = .(gene, condition)]
  return(gene_condition_means)
}

#' Reshape mean counts into wide format
#'
#' Converts a long-format table of mean gene counts into wide format,
#' with one column per condition.
#'
#' @param gene_condition_means A long-format `data.table` with columns `gene`, `condition`, and `mean_count`.
#'
#' @return A wide-format `data.table` with one row per gene and one column per condition.
#'
#' @export
reshape_means_to_wide <- function(gene_condition_means) {
  wide_gene_condition <- dcast(gene_condition_means, gene ~ condition, value.var = "mean_count")
  return(wide_gene_condition)
}

#' Bulk RNA-seq summary pipeline for wide-format counts
#'
#' Processes bulk RNA-seq counts in wide format, converts them to long format,
#' normalizes by sample totals, extracts conditions, computes per-gene means,
#' and reshapes the results to wide format for export.
#'
#' @param input_file Character string. Path to the wide-format counts file (default = `"bulk_counts_wide.csv"`).
#' @param output_file Character string. Output path for the processed summary (default = `"gene_condition_mean_counts.csv"`).
#'
#' @return Invisibly returns a wide-format `data.table` containing mean counts per gene per condition.
#'
#' @details
#' The pipeline performs the following steps:
#' \enumerate{
#'   \item Load wide-format count data.
#'   \item Convert to long format.
#'   \item Add total and normalized counts per sample.
#'   \item Extract condition labels from sample names.
#'   \item Compute per-gene mean counts by condition.
#'   \item Reshape to wide format and save the results.
#' }
#'
#' @export
task9 <- function(input_file = "bulk_counts_wide.csv",
                  output_file = "gene_condition_mean_counts.csv") {
  
  bulk_counts <- load_wide_counts(input_file)
  long_counts <- convert_wide_to_long(bulk_counts)
  long_counts <- add_sample_totals(long_counts)
  long_counts <- extract_condition(long_counts)
  
  gene_condition_means <- compute_gene_condition_means(long_counts)
  wide_gene_condition <- reshape_means_to_wide(gene_condition_means)
  
  fwrite(wide_gene_condition, output_file)
  print(head(wide_gene_condition))
  
  invisible(wide_gene_condition)
}

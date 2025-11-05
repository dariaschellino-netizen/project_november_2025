#' Identify high-expression genes and compute log2 counts
#'
#' This function reads a long-format bulk RNA count table, applies log2
#' transformation to the count values, and flags genes with high expression
#' both globally and within each gene group. It prints the first few rows of
#' the processed data for inspection and returns the modified data.table.
#'
#' @param file_path Character string. Path to the CSV file containing bulk RNA
#'   counts in long format. Default is `"bulk_counts_long.csv"`.
#' @param count_col Character string. Name of the column containing raw count
#'   values. Default is `"count"`.
#' @param group_col Character string. Name of the column defining gene groups
#'   (e.g., gene identifier). Default is `"gene"`.
#' @param global_threshold Numeric. Threshold above which counts are considered
#'   globally high. Default is `100`.
#'
#' @return A `data.table` with additional columns:
#' \describe{
#'   \item{log2_count}{Log2-transformed counts (`log2(count + 1)`).}
#'   \item{high_global}{Logical flag indicating if a count exceeds the global threshold.}
#'   \item{high_gene}{Logical flag indicating if a count is above the median within its gene group.}
#' }
#'
#' @details
#' The function uses `data.table` for efficient computation. It automatically
#' computes log2-transformed counts and two binary indicators of high expression:
#' one based on a global threshold and one relative to each gene's median.
#'
#' @examples
#' \dontrun{
#' # Using default parameters
#' dt <- task2()
#'
#' # Specifying custom file path and threshold
#' dt <- task2(file_path = "my_counts.csv", global_threshold = 150)
#' }
#'
#' @export
task2 <- function(file_path = "bulk_counts_long.csv", 
                  count_col = "count", 
                  group_col = "gene", 
                  global_threshold = 100) {
  
  dt <- fread(file_path)
  
  dt[, log2_count := log2(get(count_col) + 1)]
  
  dt[, high_global := get(count_col) > global_threshold]
  
  dt[, high_gene := get(count_col) > median(get(count_col)), by = group_col]
  
  print(head(dt))
  
  invisible(dt)
}

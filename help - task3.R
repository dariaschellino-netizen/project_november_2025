#' Read sample metadata and bulk RNA counts
#'
#' This function reads metadata and count data files into R using `data.table::fread()`
#' and returns them as a named list.
#'
#' @param metadata_file Character string. Name of the metadata CSV file located in `inst/extdata`.
#' @param counts_file Character string. Name of the counts CSV file located in `inst/extdata`.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{sample_metadata}{A `data.table` containing sample metadata.}
#'   \item{bulk_counts_long}{A `data.table` containing bulk RNA-seq counts in long format.}
#' }
#'
#' @examples
#' \dontrun{
#' files <- read_data("sample_metadata.csv", "bulk_counts_long.csv")
#' }
#'
#' @export
read_data <- function(metadata_file, counts_file) {
  metadata_path <- system.file("extdata", metadata_file, package = "library")
  counts_path <- system.file("extdata", counts_file, package = "library")

  if (!file.exists(metadata_path)) stop("Metadata file not found in package 'extdata'")
  if (!file.exists(counts_path)) stop("Counts file not found in package 'extdata'")

  sample_metadata <- data.table::fread(metadata_path)
  bulk_counts_long <- data.table::fread(counts_path)

  list(sample_metadata = sample_metadata, bulk_counts_long = bulk_counts_long)
}

#' Merge sample metadata with bulk count data
#'
#' @param bulk_counts_long A `data.table` containing bulk RNA-seq counts.
#' @param sample_metadata A `data.table` containing sample metadata.
#' @param key Character string. Column name to merge on (default: "sample_id").
#'
#' @return A merged `data.table`.
#' @export
merge_metadata <- function(bulk_counts_long, sample_metadata, key = "sample_id") {
  data.table::setkeyv(sample_metadata, key)
  merged_data <- sample_metadata[bulk_counts_long, on = key]
  return(merged_data)
}

#' Set data.table indices for efficient subsetting
#'
#' @param data A `data.table`.
#' @param indices_cols Character vector of column names to set as indices.
#'
#' @return The same `data.table` with indices set.
#' @export
set_indices <- function(data, indices_cols = c("gene", "sample_id")) {
  data.table::setindexv(data, indices_cols)
  return(data)
}

#' Benchmark subsetting performance on indexed data.tables
#'
#' @param data A `data.table` (preferably indexed).
#' @param gene_subset Character vector of genes to subset.
#' @param sample_subset Character vector of sample IDs to subset.
#' @param times Integer. Number of repetitions for benchmarking (default 10).
#'
#' @return A `microbenchmark` object.
#'
#' @examples
#' \dontrun{
#' results <- benchmark_subset(indexed_data, c("GeneA","GeneB"), c("Sample1","Sample2"))
#' print(results)
#' }
#' @export
benchmark_subset <- function(data, gene_subset, sample_subset, times = 10) {
  microbenchmark::microbenchmark(
    subset_query = data[gene %in% gene_subset & sample_id %in% sample_subset],
    times = times
  )
}

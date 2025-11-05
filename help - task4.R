# =====================================
# Task 4 Functions
# =====================================

#' Read counts and metadata for Task 4
#'
#' Loads RNA-seq counts and sample metadata files into memory using `data.table::fread()`.
#'
#' @param counts_file Character string. Path to the counts CSV file.
#' @param meta_file Character string. Path to the metadata CSV file.
#'
#' @return A named list containing:
#' \describe{
#'   \item{counts}{A `data.table` with RNA-seq count data.}
#'   \item{meta}{A `data.table` with metadata.}
#' }
#'
#' @examples
#' \dontrun{
#' data <- read_task4_data("bulk_counts_long.csv", "sample_metadata.csv")
#' }
#'
#' @export
read_task4_data <- function(counts_file, meta_file) {
  counts <- data.table::fread(counts_file)
  meta   <- data.table::fread(meta_file)
  list(counts = counts, meta = meta)
}

#' Merge metadata with count data (Task 4)
#'
#' Merges count and metadata tables by a specified key column, keeping all count records.
#'
#' @param counts A `data.table` containing RNA-seq count data.
#' @param meta A `data.table` containing sample metadata.
#' @param key Character string. Column name used for merging. Default is `"sample_id"`.
#'
#' @return A merged `data.table` containing both count and metadata information.
#'
#' @examples
#' \dontrun{
#' dt <- merge_metadata_task4(data$counts, data$meta)
#' }
#'
#' @export
merge_metadata_task4 <- function(counts, meta, key = "sample_id") {
  data.table::merge(counts, meta, by = key, all.x = TRUE)
}

#' Compute total counts per patient
#'
#' Calculates the total RNA-seq counts for each patient across all samples.
#'
#' @param dt A `data.table` containing count and metadata information.
#' @param patient_col Character string. Column name representing patient identifiers. Default is `"patient_id"`.
#' @param count_col Character string. Column name containing raw counts. Default is `"count"`.
#'
#' @return A `data.table` with one row per patient and a `total_count` column.
#'
#' @examples
#' \dontrun{
#' patient_totals <- compute_patient_totals(dt)
#' }
#'
#' @export
compute_patient_totals <- function(dt, patient_col = "patient_id", count_col = "count") {
  patient_totals <- dt[, .(total_count = sum(get(count_col))), by = get(patient_col)]
  data.table::setnames(patient_totals, "get", patient_col)  # rename grouping column
  patient_totals
}

#' Annotate data with per-patient total counts
#'
#' Merges per-patient total counts back into the full dataset.
#'
#' @param dt A `data.table` containing merged count and metadata data.
#' @param patient_totals A `data.table` returned by [compute_patient_totals()].
#' @param patient_col Character string. Column name representing patient identifiers. Default is `"patient_id"`.
#'
#' @return A `data.table` with an additional `total_count` column for each sample.
#'
#' @examples
#' \dontrun{
#' dt <- annotate_with_totals(dt, patient_totals)
#' }
#'
#' @export
annotate_with_totals <- function(dt, patient_totals, patient_col = "patient_id") {
  data.table::merge(dt, patient_totals, by = patient_col, all.x = TRUE)
}

#' Identify top expressed genes by condition
#'
#' For each experimental condition, identifies the top `n` genes
#' with the highest average expression levels.
#'
#' @param dt A `data.table` containing expression counts and metadata.
#' @param condition_col Column name for condition labels. Default is `"condition"`.
#' @param gene_col Column name for gene identifiers. Default is `"gene"`.
#' @param count_col Column name for raw counts. Default is `"count"`.
#' @param top_n Integer. Number of top genes to select per condition. Default is `10`.
#'
#' @return A `data.table` listing the top `n` genes per condition with their average counts.
#'
#' @examples
#' \dontrun{
#' top_genes <- top_genes_by_condition(dt, top_n = 10)
#' }
#'
#' @export
top_genes_by_condition <- function(dt,
                                   condition_col = "condition",
                                   gene_col = "gene",
                                   count_col = "count",
                                   top_n = 10) {
  top_genes <- dt[, .(avg_count = mean(get(count_col))), by = c(condition_col, gene_col)
  ][order(get(condition_col), -avg_count)
  ][, head(.SD, top_n), by = get(condition_col)]
  top_genes
}

#' Save Task 4 analysis results
#'
#' Writes per-patient total counts and top genes per condition to CSV files.
#'
#' @param patient_totals A `data.table` with total counts per patient.
#' @param top_genes A `data.table` with top expressed genes per condition.
#' @param patient_file Output filename for patient totals. Default `"per_patient_total_counts.csv"`.
#' @param top_genes_file Output filename for top genes. Default `"top10_genes_by_condition.csv"`.
#'
#' @return Invisibly returns `NULL` after writing files.
#'
#' @examples
#' \dontrun{
#' save_results_task4(patient_totals, top_genes)
#' }
#'
#' @export
save_results_task4 <- function(patient_totals,
                               top_genes,
                               patient_file = "per_patient_total_counts.csv",
                               top_genes_file = "top10_genes_by_condition.csv") {
  data.table::fwrite(patient_totals, patient_file)
  data.table::fwrite(top_genes, top_genes_file)
  invisible(NULL)
}

#' Task 4 main pipeline
#'
#' Runs the full Task 4 workflow: reads data, merges metadata, computes patient totals,
#' annotates samples with totals, identifies top genes per condition, and optionally saves results.
#'
#' @param counts_file Path to counts CSV file. Default `"bulk_counts_long.csv"`.
#' @param meta_file Path to metadata CSV file. Default `"sample_metadata.csv"`.
#' @param top_n Number of top genes to select per condition. Default `10`.
#' @param save_output Logical. Whether to save output CSV files. Default `TRUE`.
#'
#' @return Invisibly returns a list with:
#' \describe{
#'   \item{dt}{Merged and annotated data.table.}
#'   \item{patient_totals}{Per-patient total counts.}
#'   \item{top_genes}{Top expressed genes per condition.}
#' }
#'
#' @export
task4 <- function(counts_file = "bulk_counts_long.csv",
                  meta_file = "sample_metadata.csv",
                  top_n = 10,
                  save_output = TRUE) {
  data <- read_task4_data(counts_file, meta_file)
  dt <- merge_metadata_task4(data$counts, data$meta)
  patient_totals <- compute_patient_totals(dt)
  dt <- annotate_with_totals(dt, patient_totals)
  top_genes <- top_genes_by_condition(dt, top_n = top_n)

  if (save_output) {
    save_results_task4(patient_totals, top_genes)
  }

  print(patient_totals)
  print(top_genes)

  invisible(list(dt = dt, patient_totals = patient_totals, top_genes = top_genes))
}

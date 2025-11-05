#' Load cohort metadata files
#'
#' Loads two cohort sample metadata files (e.g., Cohort A and Cohort B) into R as `data.table` objects.
#'
#' @param fileA Character string. Path to the Cohort A samples file (default = `"cohortA_samples.csv"`).
#' @param fileB Character string. Path to the Cohort B samples file (default = `"cohortB_samples.csv"`).
#'
#' @return A list containing:
#' \describe{
#'   \item{cohortA}{A `data.table` of samples from Cohort A.}
#'   \item{cohortB}{A `data.table` of samples from Cohort B.}
#' }
#'
#' @export
load_cohorts <- function(fileA = "cohortA_samples.csv", 
                         fileB = "cohortB_samples.csv") {
  cohortA <- fread(fileA)
  cohortB <- fread(fileB)
  return(list(cohortA = cohortA, cohortB = cohortB))
}

#' Combine cohort metadata into a single dataset
#'
#' Combines metadata from two cohorts into one `data.table`, aligning columns by name.
#'
#' @param cohortA A `data.table` containing metadata for Cohort A.
#' @param cohortB A `data.table` containing metadata for Cohort B.
#'
#' @return A combined `data.table` with all samples from both cohorts.
#'
#' @export
combine_cohorts <- function(cohortA, cohortB) {
  combined <- rbindlist(list(cohortA, cohortB), use.names = TRUE, fill = TRUE)
  cat("Columns in combined dataset:\n")
  print(names(combined))
  return(combined)
}

#' Order combined cohort metadata
#'
#' Orders the combined metadata table by cohort, condition, and sample ID for clarity.
#'
#' @param combined A combined cohort `data.table` created using [combine_cohorts()].
#'
#' @return The same `data.table`, sorted by `cohort`, `condition`, and `sample_id`.
#'
#' @export
order_combined <- function(combined) {
  setorder(combined, cohort, condition, sample_id)
  return(combined)
}

#' Load bulk RNA-seq counts
#'
#' Loads a long-format RNA-seq counts table.
#'
#' @param file Character string. Path to the bulk counts file (default = `"bulk_counts_long.csv"`).
#'
#' @return A `data.table` with columns such as `gene`, `sample_id`, and `count`.
#'
#' @export
load_bulk_counts <- function(file = "bulk_counts_long.csv") {
  bulk_counts <- fread(file)
  return(bulk_counts)
}

#' Join bulk counts with cohort metadata
#'
#' Merges gene count data with cohort metadata to annotate each sample.
#'
#' @param bulk_counts A `data.table` containing long-format gene counts.
#' @param combined A `data.table` containing combined cohort metadata.
#'
#' @return A merged `data.table` including sample-level metadata and counts.
#'
#' @export
join_counts_with_metadata <- function(bulk_counts, combined) {
  counts_annot <- merge(bulk_counts, combined, by = "sample_id", all.x = TRUE)
  return(counts_annot)
}

#' Identify top variable genes
#'
#' Computes the variance of expression for each gene and selects the most variable ones.
#'
#' @param counts_annot A `data.table` containing annotated counts data with a `count` column.
#' @param top_n Integer. Number of top variable genes to select (default = 100).
#'
#' @return A character vector of gene identifiers corresponding to the top variable genes.
#'
#' @export
find_top_variable_genes <- function(counts_annot, top_n = 100) {
  gene_variances <- counts_annot[, .(variance = var(count, na.rm = TRUE)), by = gene]
  top_genes <- gene_variances[order(-variance)][1:top_n, gene]
  return(top_genes)
}

#' Filter dataset for selected top genes
#'
#' Subsets the annotated counts table to include only selected top variable genes.
#'
#' @param counts_annot A `data.table` containing annotated counts data.
#' @param top_genes A character vector of top variable genes.
#'
#' @return A filtered `data.table` containing only the selected genes.
#'
#' @export
filter_top_genes <- function(counts_annot, top_genes) {
  top_counts <- counts_annot[gene %in% top_genes]
  return(top_counts)
}

#' Compute mean gene counts by cohort and condition
#'
#' Calculates the average expression (mean count) for each gene across cohorts and conditions.
#'
#' @param top_counts A filtered `data.table` containing top genes and annotated metadata.
#'
#' @return A `data.table` summarizing mean counts for each gene–cohort–condition combination.
#'
#' @export
compute_mean_counts <- function(top_counts) {
  mean_counts <- top_counts[, .(mean_count = mean(count, na.rm = TRUE)), 
                            by = .(gene, cohort, condition)]
  return(mean_counts)
}

#' Cohort comparison pipeline
#'
#' Combines two cohorts, integrates RNA-seq counts, identifies the most variable genes,
#' and computes their mean expression per cohort and condition.
#'
#' @param fileA Path to Cohort A metadata file (default = `"cohortA_samples.csv"`).
#' @param fileB Path to Cohort B metadata file (default = `"cohortB_samples.csv"`).
#' @param bulk_counts_file Path to long-format bulk RNA-seq counts file (default = `"bulk_counts_long.csv"`).
#' @param top_n Integer specifying the number of top variable genes to analyze (default = 100).
#'
#' @return Invisibly returns a `data.table` with mean counts per gene, cohort, and condition.
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Loads two cohort metadata files.
#'   \item Combines and orders metadata.
#'   \item Loads RNA-seq counts.
#'   \item Merges counts with metadata.
#'   \item Finds the most variable genes.
#'   \item Filters and summarizes their mean counts by cohort and condition.
#' }
#' @export
task12 <- function(fileA = "cohortA_samples.csv",
                   fileB = "cohortB_samples.csv",
                   bulk_counts_file = "bulk_counts_long.csv",
                   top_n = 100) {
  
  # Load and combine cohort metadata
  cohorts <- load_cohorts(fileA, fileB)
  combined <- combine_cohorts(cohorts$cohortA, cohorts$cohortB)
  combined <- order_combined(combined)
  
  # Load and annotate bulk counts
  bulk_counts <- load_bulk_counts(bulk_counts_file)
  counts_annot <- join_counts_with_metadata(bulk_counts, combined)
  
  # Identify top variable genes and compute their mean counts
  top_genes <- find_top_variable_genes(counts_annot, top_n)
  top_counts <- filter_top_genes(counts_annot, top_genes)
  mean_counts <- compute_mean_counts(top_counts)
  
  print(head(mean_counts))
  invisible(mean_counts)
}


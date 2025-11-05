#' Load ATAC-seq and gene annotation data
#'
#' Loads ATAC-seq peak and gene annotation data from CSV files (converted from BED files).
#'
#' @param atac_file Character string. Path to the ATAC-seq peaks file (default = `"atac_peaks.bed.csv"`).
#' @param genes_file Character string. Path to the gene annotation file (default = `"gene_annotation.bed.csv"`).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{atac}{A `data.table` containing ATAC-seq peaks with columns such as `chr`, `start`, and `end`.}
#'   \item{genes}{A `data.table` containing gene annotation data with columns such as `chr`, `start`, `end`, and `gene`.}
#' }
#'
#' @export
load_task10_data <- function(atac_file = "atac_peaks.bed.csv",
                             genes_file = "gene_annotation.bed.csv") {
  atac <- fread(atac_file)      
  genes <- fread(genes_file)    
  return(list(atac = atac, genes = genes))
}

#' Ensure coordinate columns are numeric
#'
#' Converts start and end columns in ATAC and gene data to numeric types.
#'
#' @param atac A `data.table` containing ATAC peak data.
#' @param genes A `data.table` containing gene annotation data.
#'
#' @return A list containing updated `atac` and `genes` `data.table`s with numeric coordinate columns.
#'
#' @export
ensure_numeric_coords <- function(atac, genes) {
  atac[, `:=`(start = as.numeric(start), end = as.numeric(end))]
  genes[, `:=`(start = as.numeric(start), end = as.numeric(end))]
  return(list(atac = atac, genes = genes))
}

#' Set keys for overlap detection
#'
#' Sets chromosome and coordinate keys in ATAC and gene tables for use with `foverlaps()`.
#'
#' @param atac A `data.table` containing ATAC peak data with `chr`, `start`, and `end` columns.
#' @param genes A `data.table` containing gene annotation data with `chr`, `start`, and `end` columns.
#'
#' @return A list containing keyed `atac` and `genes` tables.
#'
#' @export
set_keys_for_overlap <- function(atac, genes) {
  setkey(atac, chr, start, end)
  setkey(genes, chr, start, end)
  return(list(atac = atac, genes = genes))
}

#' Find overlaps between ATAC peaks and genes
#'
#' Identifies overlapping regions between ATAC-seq peaks and gene annotations
#' using `data.table::foverlaps()`, and calculates the base pair overlap length.
#'
#' @param atac A keyed `data.table` of ATAC peaks.
#' @param genes A keyed `data.table` of gene regions.
#'
#' @return A `data.table` containing overlap pairs and an additional column `overlap_bp`
#'   representing the number of overlapping base pairs.
#'
#' @export
find_overlaps <- function(atac, genes) {
  overlaps <- foverlaps(atac, genes, nomatch = 0)
  overlaps[, overlap_bp := pmin(end, i.end) - pmax(start, i.start) + 1]
  return(overlaps)
}

#' Summarize overlap information by gene
#'
#' Aggregates overlap results to count the number of peaks and total overlap size per gene.
#'
#' @param overlaps A `data.table` produced by `find_overlaps()`, containing `gene` and `overlap_bp` columns.
#'
#' @return A `data.table` summarizing overlap statistics per gene, with columns:
#' \describe{
#'   \item{gene}{Gene name or identifier.}
#'   \item{peak_count}{Number of overlapping ATAC peaks for the gene.}
#'   \item{total_overlap_bp}{Sum of overlapping base pairs across all peaks.}
#' }
#'
#' @export
summarize_overlaps <- function(overlaps) {
  gene_overlap <- overlaps[, .(
    peak_count = .N,
    total_overlap_bp = sum(overlap_bp)
  ), by = gene]
  return(gene_overlap)
}

#' Select top genes by overlap coverage
#'
#' Returns the top genes ranked by total overlapping base pairs.
#'
#' @param gene_overlap A `data.table` containing overlap summaries per gene.
#' @param n Integer. Number of top genes to return (default = `20`).
#'
#' @return A `data.table` with the top `n` genes sorted by `total_overlap_bp` in descending order.
#'
#' @export
get_top_genes <- function(gene_overlap, n = 20) {
  top_genes <- gene_overlap[order(-total_overlap_bp)][1:n]
  return(top_genes)
}

#' ATAC–gene overlap analysis pipeline
#'
#' Runs the complete ATAC–gene overlap workflow:
#' \enumerate{
#'   \item Loads ATAC and gene data.
#'   \item Ensures numeric coordinates.
#'   \item Sets keys for overlap computation.
#'   \item Finds overlaps between peaks and genes.
#'   \item Summarizes overlaps by gene.
#'   \item Returns the top genes by total overlap length.
#' }
#'
#' @param atac_file Character string. Path to the ATAC-seq peaks file (default = `"atac_peaks.bed.csv"`).
#' @param genes_file Character string. Path to the gene annotation file (default = `"gene_annotation.bed.csv"`).
#' @param top_n Integer. Number of top genes to display (default = `20`).
#'
#' @return A `data.table` of the top `n` genes with the greatest total ATAC overlap.
#'
#' @export
task10 <- function(atac_file = "atac_peaks.bed.csv",
                   genes_file = "gene_annotation.bed.csv",
                   top_n = 20) {
  
  data_list <- load_task10_data(atac_file, genes_file)
  data_list <- ensure_numeric_coords(data_list$atac, data_list$genes)
  data_list <- set_keys_for_overlap(data_list$atac, data_list$genes)
  
  overlaps <- find_overlaps(data_list$atac, data_list$genes)
  gene_overlap <- summarize_overlaps(overlaps)
  top_genes <- get_top_genes(gene_overlap, n = top_n)
  
  print(top_genes)
}

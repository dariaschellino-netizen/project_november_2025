#' Load ATAC-seq peak data
#'
#' Reads a BED-like file (in CSV format) containing ATAC-seq peak information into a `data.table`.
#'
#' @param filepath Character string. Path to the CSV file containing ATAC-seq peaks.
#'
#' @return A `data.table` containing peak information.
#'
#' @export
load_data <- function(filepath) {
  peaks <- fread(filepath)
  return(peaks)
}

#' Filter peaks from chromosome 2 within a genomic range
#'
#' Selects ATAC-seq peaks located on a specified chromosome within given start coordinate boundaries.
#'
#' @param peaks A `data.table` containing ATAC-seq peak data.
#' @param chr_name Character string specifying the chromosome name to filter (default = `"chr2"`).
#' @param start_min Numeric. Minimum start coordinate (default = `2e6`).
#' @param start_max Numeric. Maximum start coordinate (default = `4e6`).
#'
#' @return A filtered `data.table` containing only peaks within the specified region.
#'
#' @export
filter_chr2_peaks <- function(peaks, chr_name = "chr2", start_min = 2e6, start_max = 4e6) {
  chr2_peaks <- peaks[chr == chr_name & start >= start_min & start <= start_max]
  return(chr2_peaks)
}

#' Order peaks by descending score
#'
#' Orders ATAC-seq peaks in descending order of their score (highest first).
#'
#' @param peaks A `data.table` containing ATAC-seq peaks, including a `score` column.
#'
#' @return A `data.table` sorted in descending order by score.
#'
#' @export
order_peaks_by_score <- function(peaks) {
  setorder(peaks, -score)
  return(peaks)
}

#' Select top peaks by score
#'
#' Selects the top *n* peaks with the highest scores.
#'
#' @param peaks A `data.table` of ordered peaks.
#' @param n Integer specifying how many top peaks to select (default = 50).
#'
#' @return A `data.table` containing the top *n* peaks.
#' @export
select_top_peaks <- function(peaks, n = 50) {
  top_peaks <- peaks[1:n]
  return(top_peaks)
}

#' ATAC-seq Peak Analysis Pipeline
#'
#' Loads an ATAC-seq peak dataset, filters for peaks on chromosome 2 within
#' a specific genomic range, orders them by score, and selects the top 50 peaks.
#'
#' @param filepath Character string. Path to the CSV file containing ATAC-seq peaks.
#'
#' @return A `data.table` containing the top 50 peaks from chromosome 2.
#'
#' @details
#' This function performs a simple ATAC-seq filtering workflow:
#' \enumerate{
#'   \item Loads the ATAC-seq peak file.
#'   \item Filters for peaks in the region chr2:2,000,000–4,000,000.
#'   \item Sorts by descending score.
#'   \item Selects the top 50 peaks.
#' }

#' @export
task7 <- function(filepath = "atac_peaks.bed.csv") {
  peaks <- load_data(filepath)
  chr2_peaks <- filter_chr2_peaks(peaks)
  ordered_peaks <- order_peaks_by_score(chr2_peaks)
  top50_peaks <- select_top_peaks(ordered_peaks, n = 50)
  
  print(top50_peaks)
  return(top50_peaks)
}


#' Read clinical and reference laboratory data
#'
#' Loads clinical laboratory results, reference ranges, and sample metadata into memory
#' using `data.table::fread()`.
#'
#' @param clinical_file Character string. Path to the CSV file containing clinical lab values.
#' @param reference_file Character string. Path to the CSV file containing reference ranges per lab test.
#' @param metadata_file Character string. Path to the CSV file containing patient metadata.
#'
#' @return A named list containing:
#' \describe{
#'   \item{clinical_labs}{A `data.table` with clinical lab values.}
#'   \item{lab_reference}{A `data.table` with reference ranges for each lab test.}
#'   \item{sample_metadata}{A `data.table` with sample metadata.}
#' }
#'
#'
#' @export
read_clinical_data <- function(clinical_file, reference_file, metadata_file) {
  clinical_labs <- fread(clinical_file)
  lab_reference <- fread(reference_file)
  sample_metadata <- fread(metadata_file)
  list(clinical_labs = clinical_labs, lab_reference = lab_reference, sample_metadata = sample_metadata)
}

#' Prepare clinical lab intervals for overlap classification
#'
#' Converts lab values into intervals compatible with `data.table::foverlaps()` for efficient range-based matching.
#'
#' @param clinical_labs A `data.table` containing lab test results per patient.
#' @param lab_reference A `data.table` containing reference lower and upper limits per lab test.
#'
#' @return A named list containing updated `clinical_labs` and `lab_reference` tables with interval keys set.
#'
#'
#' @export
prepare_intervals <- function(clinical_labs, lab_reference) {
  setDT(clinical_labs)
  setDT(lab_reference)
  
  clinical_labs[, `:=`(value_start = value, value_end = value)]
  
  setkey(lab_reference, lab, lower, upper)
  setkey(clinical_labs, lab, value_start, value_end)
  
  list(clinical_labs = clinical_labs, lab_reference = lab_reference)
}

#' Classify clinical lab results as normal or out-of-range
#'
#' Compares clinical lab values with reference ranges using interval overlap joins,
#' labeling each measurement as either `"normal"` or `"out_of_range"`.
#'
#' @param clinical_labs A `data.table` of lab test results with interval columns.
#' @param lab_reference A `data.table` of reference ranges per lab test.
#'
#' @return A `data.table` with all clinical lab results annotated with a `status` column.
#'
#' @details
#' The classification uses `data.table::foverlaps()` for efficient matching between
#' lab values and their reference ranges. Any lab result not overlapping with the
#' reference interval is labeled as `"out_of_range"`.
#'
#'
#' @export
classify_labs <- function(clinical_labs, lab_reference) {
  classified <- foverlaps(
    clinical_labs, lab_reference,
    by.x = c("lab", "value_start", "value_end"),
    by.y = c("lab", "lower", "upper"),
    type = "within", nomatch = 0L
  )
  classified[, status := "normal"]
  
  out_of_range <- clinical_labs[!classified, on = .(patient_id, lab, value_start, value_end)]
  out_of_range[, status := "out_of_range"]
  
  classified_labs <- rbindlist(list(classified, out_of_range), fill = TRUE)
  return(classified_labs)
}

#' Compute abnormal lab rates per patient
#'
#' Aggregates classified lab data to compute the number and proportion of abnormal lab results
#' for each patient, then merges this with sample metadata.
#'
#' @param classified_labs A `data.table` produced by [classify_labs()].
#' @param sample_metadata A `data.table` containing patient-level metadata.
#'
#' @return A `data.table` summarizing abnormal lab statistics per patient.
#'
#'
#' @export
compute_abnormal_by_patient <- function(classified_labs, sample_metadata) {
  abnormal_by_patient <- classified_labs[, .(
    total_labs = .N,
    abnormal_count = sum(status == "out_of_range"),
    abnormal_rate = sum(status == "out_of_range") / .N
  ), by = patient_id]
  
  abnormal_by_patient <- merge(abnormal_by_patient, sample_metadata, by = "patient_id", all.x = TRUE)
  return(abnormal_by_patient)
}

#' Compute abnormal lab rates per test
#'
#' Aggregates classified lab data to compute the number and proportion of abnormal
#' results per lab test across all patients.
#'
#' @param classified_labs A `data.table` produced by [classify_labs()].
#'
#' @return A `data.table` summarizing abnormal lab statistics per lab test.
#'

#' @export
compute_abnormal_by_lab <- function(classified_labs) {
  abnormal_by_lab <- classified_labs[, .(
    total_patients = .N,
    abnormal_count = sum(status == "out_of_range"),
    abnormal_rate = sum(status == "out_of_range") / .N
  ), by = lab]
  return(abnormal_by_lab)
}

#' Save classification and summary results for Task 5
#'
#' Exports the classified lab results and their per-patient and per-lab summaries as CSV files.
#'
#' @param classified_labs A `data.table` with classified lab results.
#' @param abnormal_by_patient A `data.table` with per-patient abnormal lab summaries.
#' @param abnormal_by_lab A `data.table` with per-lab abnormal lab summaries.
#' @param labs_file Character string. Output filename for classified lab results. Default is `"classified_labs.csv"`.
#' @param patient_file Character string. Output filename for per-patient results. Default is `"abnormal_by_patient.csv"`.
#' @param lab_file Character string. Output filename for per-lab results. Default is `"abnormal_by_lab.csv"`.
#'
#' @return Invisibly returns `NULL` after saving results.
#'

#'
#' @export
save_classification_results <- function(classified_labs, abnormal_by_patient, abnormal_by_lab,
                                        labs_file = "classified_labs.csv",
                                        patient_file = "abnormal_by_patient.csv",
                                        lab_file = "abnormal_by_lab.csv") {
  fwrite(classified_labs, labs_file)
  fwrite(abnormal_by_patient, patient_file)
  fwrite(abnormal_by_lab, lab_file)
}

#' Task 5 main pipeline
#'
#' Runs the complete Task 5 analysis pipeline:
#' reads clinical data and reference ranges, classifies lab results by normality,
#' computes abnormal lab rates per patient and per lab test, and optionally saves results.
#'
#' @param clinical_file Character string. Path to the clinical lab results CSV file. Default is `"clinical_labs.csv"`.
#' @param reference_file Character string. Path to the lab reference ranges CSV file. Default is `"lab_reference_ranges.csv"`.
#' @param metadata_file Character string. Path to the metadata CSV file. Default is `"sample_metadata.csv"`.
#' @param save_output Logical. Whether to save output CSV files. Default is `TRUE`.
#'
#' @return Invisibly returns a list containing:
#' \describe{
#'   \item{classified_labs}{A `data.table` with each lab result classified as normal or out-of-range.}
#'   \item{abnormal_by_patient}{Summary statistics of abnormal lab rates per patient.}
#'   \item{abnormal_by_lab}{Summary statistics of abnormal lab rates per test.}
#' }
#'
#'
#' @export
task5 <- function(clinical_file = "clinical_labs.csv",
                  reference_file = "lab_reference_ranges.csv",
                  metadata_file = "sample_metadata.csv",
                  save_output = TRUE) {
  
  data <- read_clinical_data(clinical_file, reference_file, metadata_file)
  
  prepared <- prepare_intervals(data$clinical_labs, data$lab_reference)
  
  classified_labs <- classify_labs(prepared$clinical_labs, prepared$lab_reference)
  
  abnormal_by_patient <- compute_abnormal_by_patient(classified_labs, data$sample_metadata)
  abnormal_by_lab <- compute_abnormal_by_lab(classified_labs)
  
  if (save_output) {
    save_classification_results(classified_labs, abnormal_by_patient, abnormal_by_lab)
  }
  
  print(abnormal_by_lab)
  print(abnormal_by_patient)
  
  invisible(list(classified_labs = classified_labs,
                 abnormal_by_patient = abnormal_by_patient,
                 abnormal_by_lab = abnormal_by_lab))
}


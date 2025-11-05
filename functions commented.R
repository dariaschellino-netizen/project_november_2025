#task 1
library(data.table)

# 1. Load data
load_data <- function(counts_file, metadata_file) {
  counts <- fread(counts_file)
  metadata <- fread(metadata_file)
  return(list(counts = counts, metadata = metadata))
}

# 2. Filter data
filter_data <- function(counts, metadata, gene_pattern, condition) {
  counts_filtered <- counts[gene %like% gene_pattern]
  metadata_filtered <- metadata[condition == condition]
  return(list(counts = counts_filtered, metadata = metadata_filtered))
}

# 3. Merge counts and metadata
merge_data <- function(counts_filtered, metadata_filtered) {
  DT <- merge(counts_filtered, metadata_filtered, by = "sample_id")
  return(DT)
}

# 4. Summarize data
summarize_data <- function(DT) {
  summary_dt <- DT[, .(
    mean_count = mean(count),
    median_count = median(count)
  ), by = gene]
  return(summary_dt)
}

# 5. Print top N results
print_summary <- function(summary_dt, n = 5) {
  print(head(summary_dt, n))
}

# 6. Main function to run the pipeline
task1 <- function(counts_file = "bulk_counts_long.csv",
                  metadata_file = "sample_metadata.csv",
                  gene_pattern = "^GENE_00",
                  condition = "treated") {
  data <- load_data(counts_file, metadata_file)
  filtered <- filter_data(data$counts, data$metadata, gene_pattern, condition)
  merged <- merge_data(filtered$counts, filtered$metadata)
  summary_dt <- summarize_data(merged)
  print_summary(summary_dt)
}


# Run the pipeline
task1()

#task2 

library(data.table)

task2 <- function(file_path = "bulk_counts_long.csv", 
                  count_col = "count", 
                  group_col = "gene", 
                  global_threshold = 100) {
  
  # Step 1: Load data
  dt <- fread(file_path)
  
  # Step 2: Add log2-transformed counts
  dt[, log2_count := log2(get(count_col) + 1)]
  
  # Step 3: Flag high counts using global threshold
  dt[, high_global := get(count_col) > global_threshold]
  
  # Step 4: Flag high counts using median per gene
  dt[, high_gene := get(count_col) > median(get(count_col)), by = group_col]
  
  # Step 5: Print preview
  cat("Preview of processed data:\n")
  print(head(dt))
  
  
  # Return full table invisibly
  invisible(dt)
}
#run the pipeline
task2() 

#task3
library(data.table)
library(microbenchmark)

# Function to read data
read_data <- function(metadata_file, counts_file) {
  sample_metadata <- fread(metadata_file)
  bulk_counts_long <- fread(counts_file)
  list(sample_metadata = sample_metadata, bulk_counts_long = bulk_counts_long)
}

# Function to merge metadata with counts
merge_metadata <- function(bulk_counts_long, sample_metadata, key = "sample_id") {
  setkeyv(sample_metadata, key)
  merged_data <- sample_metadata[bulk_counts_long, on = key]
  return(merged_data)
}

# Function to set indices for faster queries
set_indices <- function(data, indices_cols = c("gene", "sample_id")) {
  setindexv(data, indices_cols)
  return(data)
}

# Function to benchmark subset queries
benchmark_subset <- function(data, gene_subset, sample_subset, times = 10) {
  microbenchmark(
    subset_query = data[gene %in% gene_subset & sample_id %in% sample_subset],
    times = times
  )
}

files <- read_data("sample_metadata.csv", "bulk_counts_long.csv")

merged_data <- merge_metadata(
  bulk_counts_long = files$bulk_counts_long,
  sample_metadata = files$sample_metadata
)

indexed_data <- set_indices(merged_data)

gene_subset <- c("GeneA", "GeneB", "GeneC")
sample_subset <- c("Sample1", "Sample2")

task3 <- benchmark_subset(indexed_data, gene_subset, sample_subset)
print(task3)

#task4
library(data.table)

# --- Function to read data ---
read_task4_data <- function(counts_file, meta_file) {
  counts <- fread(counts_file)
  meta   <- fread(meta_file)
  list(counts = counts, meta = meta)
}

# --- Function to merge metadata ---
merge_metadata_task4 <- function(counts, meta, key = "sample_id") {
  dt <- merge(counts, meta, by = key, all.x = TRUE)
  return(dt)
}

# --- Function to compute per-patient total counts ---
compute_patient_totals <- function(dt, patient_col = "patient_id", count_col = "count") {
  patient_totals <- dt[, .(total_count = sum(get(count_col))), by = get(patient_col)]
  setnames(patient_totals, "get", patient_col)  # Rename the grouping column
  return(patient_totals)
}

# --- Function to annotate main table with patient totals ---
annotate_with_totals <- function(dt, patient_totals, patient_col = "patient_id") {
  dt <- merge(dt, patient_totals, by = patient_col, all.x = TRUE)
  return(dt)
}

# --- Function to get top N genes by average count per condition ---
top_genes_by_condition <- function(dt, condition_col = "condition", gene_col = "gene", count_col = "count", top_n = 10) {
  top_genes <- dt[, .(avg_count = mean(get(count_col))), by = c(condition_col, gene_col)
  ][order(get(condition_col), -avg_count)
  ][, head(.SD, top_n), by = get(condition_col)]
  return(top_genes)
}

# --- Function to save results ---
save_results_task4 <- function(patient_totals, top_genes,
                               patient_file = "per_patient_total_counts.csv",
                               top_genes_file = "top10_genes_by_condition.csv") {
  fwrite(patient_totals, patient_file)
  fwrite(top_genes, top_genes_file)
}

# --- Main pipeline function ---
task4 <- function(counts_file = "bulk_counts_long.csv", meta_file = "sample_metadata.csv",
                  top_n = 10, save_output = TRUE) {
  # Load data
  data <- read_task4_data(counts_file, meta_file)
  
  # Merge metadata
  dt <- merge_metadata_task4(data$counts, data$meta)
  
  # Compute per-patient totals
  patient_totals <- compute_patient_totals(dt)
  
  # Annotate main table with totals
  dt <- annotate_with_totals(dt, patient_totals)
  
  # Compute top genes per condition
  top_genes <- top_genes_by_condition(dt, top_n = top_n)
  
  # Save results
  if (save_output) {
    save_results_task4(patient_totals, top_genes)
  }
  
  # Print nicely
  print(patient_totals)
  print(top_genes)
  
  # Return results invisibly
  invisible(list(dt = dt, patient_totals = patient_totals, top_genes = top_genes))
}

# --- Run pipeline ---
task4()

#task5
library(data.table)

# --- Function to read data ---
read_clinical_data <- function(clinical_file, reference_file, metadata_file) {
  clinical_labs <- fread(clinical_file)
  lab_reference <- fread(reference_file)
  sample_metadata <- fread(metadata_file)
  list(clinical_labs = clinical_labs, lab_reference = lab_reference, sample_metadata = sample_metadata)
}

# --- Function to prepare intervals for foverlaps ---
prepare_intervals <- function(clinical_labs, lab_reference) {
  setDT(clinical_labs)
  setDT(lab_reference)
  
  # Create interval for each lab value
  clinical_labs[, `:=`(value_start = value, value_end = value)]
  
  # Set keys for foverlaps
  setkey(lab_reference, lab, lower, upper)
  setkey(clinical_labs, lab, value_start, value_end)
  
  list(clinical_labs = clinical_labs, lab_reference = lab_reference)
}

# --- Function to classify lab values ---
classify_labs <- function(clinical_labs, lab_reference) {
  # Interval join
  classified <- foverlaps(
    clinical_labs, lab_reference,
    by.x = c("lab", "value_start", "value_end"),
    by.y = c("lab", "lower", "upper"),
    type = "within", nomatch = 0L
  )
  classified[, status := "normal"]
  
  # Identify out-of-range values
  out_of_range <- clinical_labs[!classified, on = .(patient_id, lab, value_start, value_end)]
  out_of_range[, status := "out_of_range"]
  
  # Combine both
  classified_labs <- rbindlist(list(classified, out_of_range), fill = TRUE)
  return(classified_labs)
}

# --- Function to compute abnormal rates by patient ---
compute_abnormal_by_patient <- function(classified_labs, sample_metadata) {
  abnormal_by_patient <- classified_labs[, .(
    total_labs = .N,
    abnormal_count = sum(status == "out_of_range"),
    abnormal_rate = sum(status == "out_of_range") / .N
  ), by = patient_id]
  
  # Merge metadata
  abnormal_by_patient <- merge(abnormal_by_patient, sample_metadata, by = "patient_id", all.x = TRUE)
  return(abnormal_by_patient)
}

# --- Function to compute abnormal rates by lab ---
compute_abnormal_by_lab <- function(classified_labs) {
  abnormal_by_lab <- classified_labs[, .(
    total_patients = .N,
    abnormal_count = sum(status == "out_of_range"),
    abnormal_rate = sum(status == "out_of_range") / .N
  ), by = lab]
  return(abnormal_by_lab)
}

# --- Function to save results ---
save_classification_results <- function(classified_labs, abnormal_by_patient, abnormal_by_lab,
                                        labs_file = "classified_labs.csv",
                                        patient_file = "abnormal_by_patient.csv",
                                        lab_file = "abnormal_by_lab.csv") {
  fwrite(classified_labs, labs_file)
  fwrite(abnormal_by_patient, patient_file)
  fwrite(abnormal_by_lab, lab_file)
}

# --- Main pipeline function ---
task5 <- function(clinical_file = "clinical_labs.csv",
                                  reference_file = "lab_reference_ranges.csv",
                                  metadata_file = "sample_metadata.csv",
                                  save_output = TRUE) {
  
  # Load data
  data <- read_clinical_data(clinical_file, reference_file, metadata_file)
  
  # Prepare intervals
  prepared <- prepare_intervals(data$clinical_labs, data$lab_reference)
  
  # Classify labs
  classified_labs <- classify_labs(prepared$clinical_labs, prepared$lab_reference)
  
  # Compute abnormal rates
  abnormal_by_patient <- compute_abnormal_by_patient(classified_labs, data$sample_metadata)
  abnormal_by_lab <- compute_abnormal_by_lab(classified_labs)
  
  # Save results
  if (save_output) {
    save_classification_results(classified_labs, abnormal_by_patient, abnormal_by_lab)
  }
  
  # Print summary
  print(abnormal_by_lab)
  print(abnormal_by_patient)
  
  invisible(list(classified_labs = classified_labs,
                 abnormal_by_patient = abnormal_by_patient,
                 abnormal_by_lab = abnormal_by_lab))
}

# --- Run pipeline ---
task5()

#task7

library(data.table)

# --- Function 1: Load data ---
load_data <- function(filepath) {
  peaks <- fread(filepath)
  return(peaks)
}

# --- Function 2: Filter peaks for chr2 and coordinate range ---
filter_chr2_peaks <- function(peaks, chr_name = "chr2", start_min = 2e6, start_max = 4e6) {
  # Compare the 'chr' column to chr_name string
  chr2_peaks <- peaks[chr == chr_name & start >= start_min & start <= start_max]
  return(chr2_peaks)
}

# --- Function 3: Order peaks by score ---
order_peaks_by_score <- function(peaks) {
  setorder(peaks, -score)
  return(peaks)
}

# --- Function 4: Select top N peaks ---
select_top_peaks <- function(peaks, n = 50) {
  top_peaks <- peaks[1:n]
  return(top_peaks)
}

# --- Main pipeline function ---
task7 <- function(filepath = "atac_peaks.bed.csv") {
  peaks <- load_data(filepath)
  chr2_peaks <- filter_chr2_peaks(peaks)
  ordered_peaks <- order_peaks_by_score(chr2_peaks)
  top50_peaks <- select_top_peaks(ordered_peaks, n = 50)
  
  print(top50_peaks)
  return(top50_peaks)
}

# run pipeline
# task7()

#task8
library(data.table)

# --- Function 1: Load data ---
load_task8_data <- function(counts_file = "bulk_counts_long.csv",
                            metadata_file = "sample_metadata.csv") {
  counts <- fread(counts_file)        # Expect columns: gene, sample_id, count
  metadata <- fread(metadata_file)    # Expect columns: sample_id, condition
  return(list(counts = counts, metadata = metadata))
}

# --- Function 2: Merge counts with metadata ---
merge_counts_metadata <- function(counts, metadata) {
  dt <- merge(counts, metadata, by = "sample_id")
  return(dt)
}

# --- Function 3: Compute per-condition summary statistics per gene ---
compute_summary_stats <- function(dt) {
  summary_dt <- dt[, .(
    mean_count = mean(count, na.rm = TRUE),
    median_count = median(count, na.rm = TRUE),
    Q1 = quantile(count, 0.25, na.rm = TRUE),
    Q3 = quantile(count, 0.75, na.rm = TRUE)
  ), by = .(gene, condition)]
  return(summary_dt)
}

# --- Function 4: Reshape to wide format ---
reshape_summary_wide <- function(summary_dt, value_col = "mean_count") {
  summary_wide <- dcast(summary_dt, gene ~ condition, value.var = value_col)
  return(summary_wide)
}

# --- Function 5: Filter genes based on fold-change ---
filter_genes_fc <- function(summary_wide, control_name = "control", treated_name = "treated", fold_change = 2) {
  filtered_genes <- summary_wide[get(treated_name) >= fold_change * get(control_name),
                                 .(gene, control = get(control_name), treated = get(treated_name))]
  return(filtered_genes)
}

# --- Function 6: Merge other stats back for filtered genes ---
merge_other_stats <- function(filtered_genes, summary_dt) {
  other_stats <- dcast(summary_dt, gene ~ condition, value.var = c("median_count", "Q1", "Q3"))
  result <- merge(filtered_genes, other_stats, by = "gene")
  return(result)
}

# --- Main pipeline function ---
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
  
}

# run the pipeline
# task8()


#task9
library(data.table)

# --- Function 1: Load wide-count data ---
load_wide_counts <- function(filepath = "bulk_counts_wide.csv") {
  bulk_counts <- fread(filepath)
  return(bulk_counts)
}

# --- Function 2: Convert wide to long format ---
convert_wide_to_long <- function(bulk_counts, id_col = "gene", value_name = "count") {
  long_counts <- melt(
    bulk_counts,
    id.vars = id_col,
    variable.name = "sample",
    value.name = value_name
  )
  return(long_counts)
}

# --- Function 3: Add per-sample totals and normalized counts ---
add_sample_totals <- function(long_counts) {
  sample_totals <- long_counts[, .(total_count = sum(count)), by = sample]
  long_counts <- merge(long_counts, sample_totals, by = "sample")
  long_counts[, normalized_count := count / total_count]
  return(long_counts)
}

# --- Function 4: Extract condition from sample names ---
extract_condition <- function(long_counts, pattern_split = "_", part = 1) {
  long_counts[, condition := tstrsplit(sample, pattern_split)[[part]]]
  return(long_counts)
}

# --- Function 5: Compute mean counts per gene per condition ---
compute_gene_condition_means <- function(long_counts) {
  gene_condition_means <- long_counts[, .(mean_count = mean(count, na.rm = TRUE)),
                                      by = .(gene, condition)]
  return(gene_condition_means)
}

# --- Function 6: Convert long summary to wide format ---
reshape_means_to_wide <- function(gene_condition_means) {
  wide_gene_condition <- dcast(gene_condition_means, gene ~ condition, value.var = "mean_count")
  return(wide_gene_condition)
}

# --- Main pipeline function ---
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
  
}

# run the pipeline
# task9()

#task10

library(data.table)

# --- Function 1: Load data ---
load_task10_data <- function(atac_file = "atac_peaks.bed.csv",
                             genes_file = "gene_annotation.bed.csv") {
  atac <- fread(atac_file)      # columns: chr, start, end
  genes <- fread(genes_file)    # columns: chr, start, end, gene_id
  return(list(atac = atac, genes = genes))
}

# --- Function 2: Ensure numeric coordinates ---
ensure_numeric_coords <- function(atac, genes) {
  atac[, `:=`(start = as.numeric(start), end = as.numeric(end))]
  genes[, `:=`(start = as.numeric(start), end = as.numeric(end))]
  return(list(atac = atac, genes = genes))
}

# --- Function 3: Set keys for overlap join ---
set_keys_for_overlap <- function(atac, genes) {
  setkey(atac, chr, start, end)
  setkey(genes, chr, start, end)
  return(list(atac = atac, genes = genes))
}

# --- Function 4: Find overlaps using foverlaps ---
find_overlaps <- function(atac, genes) {
  overlaps <- foverlaps(atac, genes, nomatch = 0)
  overlaps[, overlap_bp := pmin(end, i.end) - pmax(start, i.start) + 1]
  return(overlaps)
}

# --- Function 5: Summarize overlaps per gene ---
summarize_overlaps <- function(overlaps) {
  gene_overlap <- overlaps[, .(
    peak_count = .N,
    total_overlap_bp = sum(overlap_bp)
  ), by = gene]
  return(gene_overlap)
}

# --- Function 6: Get top N genes by total overlap ---
get_top_genes <- function(gene_overlap, n = 20) {
  top_genes <- gene_overlap[order(-total_overlap_bp)][1:n]
  return(top_genes)
}

# --- Main pipeline function ---
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
  return(top_genes)
}

# run the pipeline
# task10()

#task11
library(data.table)

# --- Function 1: Load input files ---
load_task11_data <- function(variants_file = "variants.csv",
                             genes_file = "gene_annotation.bed.csv") {
  variants <- fread(variants_file)   # columns: sample_id, chr, pos, impact
  genes <- fread(genes_file)         # columns: chr, start, end, gene_name
  return(list(variants = variants, genes = genes))
}

# --- Function 2: Convert SNPs to 1-bp intervals ---
convert_snps_to_intervals <- function(variants) {
  variants[, start := pos]
  variants[, end := pos]
  return(variants)
}

# --- Function 3: Set keys for foverlaps ---
set_keys_for_overlap <- function(variants, genes) {
  setkey(variants, chr, start, end)
  setkey(genes, chr, start, end)
  return(list(variants = variants, genes = genes))
}

# --- Function 4: Find overlaps between variants and genes ---
find_variant_gene_overlaps <- function(variants, genes) {
  variant_gene_map <- foverlaps(variants, genes, nomatch = 0L)
  return(variant_gene_map)
}

# --- Function 5: Filter HIGH-impact variants ---
filter_high_impact <- function(variant_gene_map) {
  high_impact <- variant_gene_map[impact == "HIGH"]
  return(high_impact)
}

# --- Function 6: Summarize HIGH-impact counts per gene and sample ---
summarize_high_impact_counts <- function(high_impact) {
  summary_counts <- high_impact[, .N, by = .(gene_name, sample_id)]
  setnames(summary_counts, "N", "HIGH_impact_count")
  return(summary_counts)
}

# --- Function 7: List genes with HIGH-impact variants ---
list_genes_with_high_impact <- function(high_impact) {
  genes_with_high_impact <- unique(high_impact$gene_name)
  genes_with_high_impact_dt <- data.table(gene_name = genes_with_high_impact)
  return(genes_with_high_impact_dt)
}

# --- Main pipeline function ---
task11 <- function(variants_file = "variants.csv",
                   genes_file = "gene_annotation.bed.csv") {
  
  # Step 1: Load data
  data_list <- load_task11_data(variants_file, genes_file)
  
  # Step 2: Convert SNPs to 1-bp intervals
  variants <- convert_snps_to_intervals(data_list$variants)
  
  # Step 3: Set keys
  data_list <- set_keys_for_overlap(variants, data_list$genes)
  
  # Step 4: Find overlaps
  variant_gene_map <- find_variant_gene_overlaps(data_list$variants, data_list$genes)
  
  # Step 5: Filter HIGH-impact variants
  high_impact <- filter_high_impact(variant_gene_map)
  
  # Step 6: Summarize counts per gene and sample
  summary_counts <- summarize_high_impact_counts(high_impact)
  
  # Step 7: List genes with HIGH-impact variants
  genes_with_high_impact_dt <- list_genes_with_high_impact(high_impact)
  
  # Step 8: Save results
  fwrite(summary_counts, "high_impact_summary.csv")
  fwrite(genes_with_high_impact_dt, "genes_with_high_impact.csv")
  
  # Quick check
  print(summary_counts)
  print(genes_with_high_impact_dt)
  
  
  
}

# run the pipeline
# task11()


#task12
library(data.table)

# --- Function 1: Load cohort data ---
load_cohorts <- function(fileA = "cohortA_samples.csv", 
                         fileB = "cohortB_samples.csv") {
  cohortA <- fread(fileA)
  cohortB <- fread(fileB)
  return(list(cohortA = cohortA, cohortB = cohortB))
}

# --- Function 2: Combine cohorts safely ---
combine_cohorts <- function(cohortA, cohortB) {
  combined <- rbindlist(list(cohortA, cohortB), use.names = TRUE, fill = TRUE)
  cat("Columns in combined dataset:\n")
  print(names(combined))
  return(combined)
}

# --- Function 3: Order combined data ---
order_combined <- function(combined) {
  setorder(combined, cohort, condition, sample_id)
  return(combined)
}

# --- Function 4: Load bulk counts in long format ---
load_bulk_counts <- function(file = "bulk_counts_long.csv") {
  bulk_counts <- fread(file)
  return(bulk_counts)
}

# --- Function 5: Join combined metadata with counts ---
join_counts_with_metadata <- function(bulk_counts, combined) {
  counts_annot <- merge(bulk_counts, combined, by = "sample_id", all.x = TRUE)
  return(counts_annot)
}

# --- Function 6: Identify top N most variable genes ---
find_top_variable_genes <- function(counts_annot, top_n = 100) {
  gene_variances <- counts_annot[, .(variance = var(count, na.rm = TRUE)), by = gene]
  top_genes <- gene_variances[order(-variance)][1:top_n, gene]
  return(top_genes)
}

# --- Function 7: Filter counts for top genes ---
filter_top_genes <- function(counts_annot, top_genes) {
  top_counts <- counts_annot[gene %in% top_genes]
  return(top_counts)
}

# --- Function 8: Compute per-cohort, per-condition mean counts ---
compute_mean_counts <- function(top_counts) {
  mean_counts <- top_counts[, .(mean_count = mean(count, na.rm = TRUE)), 
                            by = .(gene, cohort, condition)]
  return(mean_counts)
}

# --- Main pipeline function ---
task12 <- function(fileA = "cohortA_samples.csv",
                   fileB = "cohortB_samples.csv",
                   bulk_counts_file = "bulk_counts_long.csv",
                   top_n = 100) {
  
  # Step 1–3: Load and combine cohorts
  cohorts <- load_cohorts(fileA, fileB)
  combined <- combine_cohorts(cohorts$cohortA, cohorts$cohortB)
  combined <- order_combined(combined)
  
  # Step 4–5: Load counts and join metadata
  bulk_counts <- load_bulk_counts(bulk_counts_file)
  counts_annot <- join_counts_with_metadata(bulk_counts, combined)
  
  # Step 6–8: Find top variable genes and compute summary
  top_genes <- find_top_variable_genes(counts_annot, top_n)
  top_counts <- filter_top_genes(counts_annot, top_genes)
  mean_counts <- compute_mean_counts(top_counts)
  
  cat("\nMean counts per gene per cohort per condition (first 6 rows):\n")
  print(head(mean_counts))
  
}

# run the pipeline
# task12()

#final revision
#-----------------------------------
# Final Revision: Cell Type & Integration Cluster Analysis (modular version)
#-----------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)

# --- Function 1: Load data ---
load_final_revision_data <- function(integration_file = "annotated_GSM3516673_normal_annotated_GSM3516672_tumor_SeuratIntegration.csv",
                                     celltype_file = "nt_combined_clustering.output.csv") {
  integration_df <- read_csv(integration_file, show_col_types = FALSE)
  celltype_df <- read_csv(celltype_file, show_col_types = FALSE)
  return(list(integration_df = integration_df, celltype_df = celltype_df))
}

# --- Function 2: Clean and align cell IDs ---
clean_cell_ids <- function(integration_df, celltype_df) {
  integration_df$cell <- gsub("_X_", "", integration_df$cell)
  integration_df$cell <- str_trim(tolower(integration_df$cell))
  celltype_df$cell <- str_trim(tolower(celltype_df$cell))
  return(list(integration_df = integration_df, celltype_df = celltype_df))
}

# --- Function 3: Merge datasets ---
merge_datasets <- function(integration_df, celltype_df) {
  merged_df <- inner_join(celltype_df, integration_df, by = "cell")
  message(paste("Merged dataset rows:", nrow(merged_df)))
  write_csv(merged_df, "merged_celltype_clusters.csv")
  return(merged_df)
}

# --- Function 4: Count cells per cell type per cluster ---
count_cells_per_cluster <- function(merged_df) {
  count_df <- merged_df %>%
    group_by(integration_cluster, cell_type) %>%
    summarise(count = n(), .groups = "drop")
  write_csv(count_df, "celltype_per_cluster_counts.csv")
  return(count_df)
}

# --- Function 5: Create summary table with tissue info ---
create_summary_table <- function(merged_df) {
  summary_df <- merged_df %>%
    group_by(integration_cluster, cell_type, sample_type) %>%
    summarise(count = n(), .groups = "drop")
  write_csv(summary_df, "celltype_cluster_tissue_summary.csv")
  return(summary_df)
}

# --- Function 6: Prepare safe plotting data ---
prepare_plot_data <- function(merged_df) {
  merged_df$sample_type <- factor(merged_df$sample_type, levels = c("N", "T"))
  plot_df <- merged_df %>%
    group_by(integration_cluster, sample_type, cell_type) %>%
    summarise(count = n(), .groups = "drop") %>%
    tidyr::complete(integration_cluster, sample_type, cell_type, fill = list(count = 0))
  return(plot_df)
}

# --- Function 7: Plot distribution of cell types per cluster by tissue ---
plot_distribution <- function(plot_df) {
  ggplot(plot_df, aes(x = integration_cluster, y = count, fill = cell_type)) +
    geom_bar(stat = "identity") +
    facet_wrap(~sample_type, scales = "free_y") +
    labs(
      title = "Cell Type Distribution per Cluster by Tissue",
      x = "Integration Cluster",
      y = "Number of Cells"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set3")
}

# --- Function 8: Calculate normalized percentages per cluster per tissue ---
calculate_percentages <- function(merged_df) {
  percent_df <- merged_df %>%
    group_by(integration_cluster, sample_type, cell_type) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(integration_cluster, sample_type) %>%
    mutate(percent = count / sum(count) * 100) %>%
    ungroup() %>%
    tidyr::complete(integration_cluster, sample_type, cell_type, fill = list(percent = 0))
  
  write_csv(percent_df, "celltype_cluster_tissue_percent.csv")
  return(percent_df)
}

# --- Function 9: Plot normalized percentages ---
plot_percentages <- function(percent_df) {
  ggplot(percent_df, aes(x = integration_cluster, y = percent, fill = cell_type)) +
    geom_bar(stat = "identity") +
    facet_wrap(~sample_type, scales = "free_y") +
    labs(
      title = "Normalized % of Cell Types per Cluster by Tissue",
      x = "Integration Cluster",
      y = "Percentage of Cells"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set3")
}

# --- Main pipeline function ---
final_revision <- function(integration_file = "annotated_GSM3516673_normal_annotated_GSM3516672_tumor_SeuratIntegration.csv",
                           celltype_file = "nt_combined_clustering.output.csv") {
  
  # Load and clean data
  data_list <- load_final_revision_data(integration_file, celltype_file)
  data_list <- clean_cell_ids(data_list$integration_df, data_list$celltype_df)
  
  # Merge datasets
  merged_df <- merge_datasets(data_list$integration_df, data_list$celltype_df)
  
  # Summaries
  count_df <- count_cells_per_cluster(merged_df)
  summary_df <- create_summary_table(merged_df)
  
  # Plot 1: Raw counts per cluster
  plot_df <- prepare_plot_data(merged_df)
  print(plot_distribution(plot_df))
  
  # Plot 2: Percent normalized per cluster
  percent_df <- calculate_percentages(merged_df)
  print(plot_percentages(percent_df))
  
  message("Final revision pipeline completed successfully.")
  
  
}

# run the pipeline
# final_revision()



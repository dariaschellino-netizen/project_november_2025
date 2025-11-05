#' Load integration and cell type annotation data
#'
#' Loads two CSV files: one containing integrated Seurat data and one containing
#' cell type clustering results.
#'
#' @param integration_file Character string. Path to the Seurat integration results file (default: `"annotated_GSM3516673_normal_annotated_GSM3516672_tumor_SeuratIntegration.csv"`).
#' @param celltype_file Character string. Path to the cell type clustering results file (default: `"nt_combined_clustering.output.csv"`).
#'
#' @return A list with:
#' \describe{
#'   \item{integration_df}{A tibble with integrated Seurat data.}
#'   \item{celltype_df}{A tibble with cell type cluster annotations.}
#' }
#'
#' @export
load_final_revision_data <- function(integration_file = "annotated_GSM3516673_normal_annotated_GSM3516672_tumor_SeuratIntegration.csv",
                                     celltype_file = "nt_combined_clustering.output.csv") {
  integration_df <- read_csv(integration_file, show_col_types = FALSE)
  celltype_df <- read_csv(celltype_file, show_col_types = FALSE)
  return(list(integration_df = integration_df, celltype_df = celltype_df))
}

#' Clean cell IDs for consistency between datasets
#'
#' Standardizes cell IDs by removing unwanted substrings and ensuring consistent casing.
#'
#' @param integration_df A tibble containing Seurat integration results.
#' @param celltype_df A tibble containing cell type clustering data.
#'
#' @return A list with cleaned versions of both input data frames.
#'
#' @export
clean_cell_ids <- function(integration_df, celltype_df) {
  integration_df$cell <- gsub("_X_", "", integration_df$cell)
  integration_df$cell <- str_trim(tolower(integration_df$cell))
  celltype_df$cell <- str_trim(tolower(celltype_df$cell))
  return(list(integration_df = integration_df, celltype_df = celltype_df))
}

#' Merge integration and cell type datasets
#'
#' Merges the Seurat integration results with cell type clustering annotations
#' using an inner join on the `cell` column.
#'
#' @param integration_df A tibble of integrated Seurat data.
#' @param celltype_df A tibble of cell type annotations.
#'
#' @return A merged tibble containing both integration and cell type information.
#'
#' @export
merge_datasets <- function(integration_df, celltype_df) {
  merged_df <- inner_join(celltype_df, integration_df, by = "cell")
  message(paste("Merged dataset rows:", nrow(merged_df)))
  write_csv(merged_df, "merged_celltype_clusters.csv")
  return(merged_df)
}

#' Count cells per integration cluster and cell type
#'
#' Counts how many cells of each type are present in each integration cluster.
#'
#' @param merged_df A merged tibble returned by [merge_datasets()].
#'
#' @return A tibble with columns: `integration_cluster`, `cell_type`, and `count`.
#'
#' @export
count_cells_per_cluster <- function(merged_df) {
  count_df <- merged_df %>%
    group_by(integration_cluster, cell_type) %>%
    summarise(count = n(), .groups = "drop")
  
  write_csv(count_df, "celltype_per_cluster_counts.csv")
  return(count_df)
}


#' Create a summary table including sample type
#'
#' Summarizes the number of cells per cluster, cell type, and tissue sample type.
#'
#' @param merged_df A merged tibble containing cluster, cell type, and sample type columns.
#'
#' @return A tibble summarizing counts by cluster, cell type, and sample type.
#'
#' @export
create_summary_table <- function(merged_df) {
  summary_df <- merged_df %>%
    group_by(integration_cluster, cell_type, sample_type) %>%
    summarise(count = n(), .groups = "drop")
  write_csv(summary_df, "celltype_cluster_tissue_summary.csv")
  return(summary_df)
}

#' Prepare data for visualization
#'
#' Prepares count data for plotting by ensuring all combinations of cluster, sample type,
#' and cell type are present (missing combinations are filled with zero counts).
#'
#' @param merged_df A merged tibble containing `integration_cluster`, `sample_type`, and `cell_type`.
#'
#' @return A complete tibble suitable for plotting.
#'
#' @export
prepare_plot_data <- function(merged_df) {
  merged_df$sample_type <- factor(merged_df$sample_type, levels = c("N", "T"))
  plot_df <- merged_df %>%
    group_by(integration_cluster, sample_type, cell_type) %>%
    summarise(count = n(), .groups = "drop") %>%
    tidyr::complete(integration_cluster, sample_type, cell_type, fill = list(count = 0))
  return(plot_df)
}

#' Plot raw counts distribution
#'
#' Creates a stacked bar plot showing raw cell counts per integration cluster,
#' faceted by sample type.
#'
#' @param plot_df A tibble created by [prepare_plot_data()].
#'
#' @return A `ggplot` object visualizing cell type distributions.
#'
#' @export
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

#' Calculate percentage-based distributions
#'
#' Computes normalized percentages of cell types within each integration cluster
#' and sample type.
#'
#' @param merged_df A merged tibble containing integration cluster, sample type, and cell type.
#'
#' @return A tibble with columns: `integration_cluster`, `sample_type`, `cell_type`, and `percent`.
#'
#' @export
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

#' Plot percentage-based distributions
#'
#' Creates a stacked bar plot showing normalized cell-type proportions per integration cluster.
#'
#' @param percent_df A tibble produced by [calculate_percentages()].
#'
#' @return A `ggplot` object showing relative percentages per cluster and tissue type.
#'
#' @export
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

#' Final analysis and visualization pipeline
#'
#' A complete workflow to load Seurat integration data and cell type annotations,
#' clean and merge datasets, compute counts and percentages, and generate publication-ready plots.
#'
#' @param integration_file Path to the integrated Seurat CSV file.
#' @param celltype_file Path to the cell type clustering CSV file.
#'
#' @return Invisibly returns a list of generated data frames and plots.
#'
#' @details
#' Steps performed:
#' \enumerate{
#'   \item Load integration and cell type data.
#'   \item Clean and standardize cell IDs.
#'   \item Merge datasets.
#'   \item Summarize counts per cluster and cell type.
#'   \item Compute percentages for normalization.
#'   \item Generate two plots: raw counts and percentages.
#' }
#'
#' @export
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
  
  invisible(list(
    merged_df = merged_df,
    count_df = count_df,
    summary_df = summary_df,
    plot_df = plot_df,
    percent_df = percent_df
  ))
}
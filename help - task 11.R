#' Load variant and gene annotation data
#'
#' Loads variant (e.g., SNP or indel) and gene annotation data from CSV files.
#'
#' @param variants_file Character string. Path to the variants file (default = `"variants.csv"`).
#' @param genes_file Character string. Path to the gene annotation file (default = `"gene_annotation.bed.csv"`).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{variants}{A `data.table` containing variant data with columns such as `chr`, `pos`, `sample_id`, and `impact`.}
#'   \item{genes}{A `data.table` containing gene annotation data with columns such as `chr`, `start`, `end`, and `gene`.}
#' }
#'
#' @export
load_task11_data <- function(variants_file = "variants.csv",
                             genes_file = "gene_annotation.bed.csv") {
  variants <- fread(variants_file)  
  genes <- fread(genes_file)         
  return(list(variants = variants, genes = genes))
}

#' Convert SNP positions to interval format
#'
#' Converts variant positions (single base coordinates) into interval format
#' compatible with `data.table::foverlaps()` for overlap operations.
#'
#' @param variants A `data.table` containing variant data with a `pos` column.
#'
#' @return A `data.table` identical to the input but with additional columns `start` and `end`
#'   set equal to `pos`.
#'
#' @export
convert_snps_to_intervals <- function(variants) {
  variants[, start := pos]
  variants[, end := pos]
  return(variants)
}

#' Set keys for overlap computation
#'
#' Sets chromosome and coordinate keys for `variants` and `genes` tables to enable overlap joins.
#'
#' @param variants A `data.table` containing variants with columns `chr`, `start`, and `end`.
#' @param genes A `data.table` containing gene regions with columns `chr`, `start`, and `end`.
#'
#' @return A list with two keyed `data.table` objects: `variants` and `genes`.
#'
#' @export
set_keys_for_overlap <- function(variants, genes) {
  setkey(variants, chr, start, end)
  setkey(genes, chr, start, end)
  return(list(variants = variants, genes = genes))
}

#' Find overlaps between variants and genes
#'
#' Identifies overlaps between genomic variants and gene regions using
#' `data.table::foverlaps()`.
#'
#' @param variants A keyed `data.table` of variants.
#' @param genes A keyed `data.table` of genes.
#'
#' @return A `data.table` where each variant is annotated with the gene(s)
#' it overlaps.
#'
#' @export
find_variant_gene_overlaps <- function(variants, genes) {
  variant_gene_map <- foverlaps(variants, genes, nomatch = 0L)
  return(variant_gene_map)
}

#' Filter variants with high predicted impact
#'
#' Filters overlapping variants to retain only those with `"HIGH"` impact annotations.
#'
#' @param variant_gene_map A `data.table` containing variant–gene overlap results with an `impact` column.
#'
#' @return A subset `data.table` containing only high-impact variants.
#'
#' @export
filter_high_impact <- function(variant_gene_map) {
  high_impact <- variant_gene_map[impact == "HIGH"]
  return(high_impact)
}

#' Summarize high-impact variant counts
#'
#' Counts the number of high-impact variants per gene and per sample.
#'
#' @param high_impact A `data.table` containing high-impact variants with `gene` and `sample_id` columns.
#'
#' @return A `data.table` with columns:
#' \describe{
#'   \item{gene}{Gene identifier.}
#'   \item{sample_id}{Sample identifier.}
#'   \item{HIGH_impact_count}{Number of high-impact variants for the gene in that sample.}
#' }
#'
#' @export
summarize_high_impact_counts <- function(high_impact) {
  summary_counts <- high_impact[, .N, by = .(gene, sample_id)]
  setnames(summary_counts, "N", "HIGH_impact_count")
  return(summary_counts)
}

#' List genes containing high-impact variants
#'
#' Extracts a list of unique genes containing at least one high-impact variant.
#'
#' @param high_impact A `data.table` containing high-impact variant annotations with a `gene` column.
#'
#' @return A `data.table` with a single column:
#' \describe{
#'   \item{gene}{Unique gene names with at least one high-impact variant.}
#' }
#'
#' @export
list_genes_with_high_impact <- function(high_impact) {
  genes_with_high_impact <- unique(high_impact$gene)
  genes_with_high_impact_dt <- data.table(gene = genes_with_high_impact)
  return(genes_with_high_impact_dt)
}

#' Variant–gene overlap and impact summary pipeline
#'
#' Runs the complete workflow to identify high-impact variants and summarize their
#' distribution across genes and samples.
#'
#' The steps include:
#' \enumerate{
#'   \item Loading variant and gene data.
#'   \item Converting SNPs to interval format.
#'   \item Setting overlap keys.
#'   \item Finding variant–gene overlaps.
#'   \item Filtering high-impact variants.
#'   \item Summarizing counts and listing affected genes.
#' }
#'
#' @param variants_file Character string. Path to the variants file (default = `"variants.csv"`).
#' @param genes_file Character string. Path to the gene annotation file (default = `"gene_annotation.bed.csv"`).
#'
#' @return Invisibly returns a list containing:
#' \describe{
#'   \item{summary_counts}{A summary table of high-impact variant counts per gene and sample.}
#'   \item{genes_with_high_impact_dt}{A table listing genes with at least one high-impact variant.}
#' }
#'
#' @details
#' Two CSV files are also written to disk:
#' \itemize{
#'   \item `"high_impact_summary.csv"` — summary of high-impact variants per gene and sample.
#'   \item `"genes_with_high_impact.csv"` — list of genes with at least one high-impact variant.
#' }
#' @export
task11 <- function(variants_file = "variants.csv",
                   genes_file = "gene_annotation.bed.csv") {
  
  # Load data
  data_list <- load_task11_data(variants_file, genes_file)
  
  # Convert SNP positions to intervals
  variants <- convert_snps_to_intervals(data_list$variants)
  
  # Set keys for overlap
  data_list <- set_keys_for_overlap(variants, data_list$genes)
  
  # Find overlaps
  variant_gene_map <- find_variant_gene_overlaps(data_list$variants, data_list$genes)
  
  # Filter for high-impact variants
  high_impact <- filter_high_impact(variant_gene_map)
  
  # Summarize counts
  summary_counts <- summarize_high_impact_counts(high_impact)
  
  # List genes with high-impact variants
  genes_with_high_impact_dt <- list_genes_with_high_impact(high_impact)
  
  # Save outputs
  fwrite(summary_counts, "high_impact_summary.csv")
  fwrite(genes_with_high_impact_dt, "genes_with_high_impact.csv")
  
  # Display summaries
  print(summary_counts)
  print(genes_with_high_impact_dt)
  
  invisible(list(
    summary_counts = summary_counts,
    genes_with_high_impact_dt = genes_with_high_impact_dt
  ))
}

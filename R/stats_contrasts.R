#' Write statistical results tables to disk
#'
#' Writes `all_results.tsv`, `significant_results.tsv`, `upregulated.tsv`,
#' and `downregulated.tsv` for a single contrast.
#'
#' @param result_dt  `data.table` from `run_limma_deqms()` for one contrast.
#' @param output_dir Directory to write TSV files.
#' @param project    Project name for filenames.
#' @param contrast   Contrast name string.
#' @param level      Analysis level string (e.g. `"protein"`).
#'
#' @return Named character vector of written file paths (invisibly).
#'
#' @examples
#' \dontrun{
#' write_contrast_tables(result_dt, "results/04_statistics/tables",
#'                       "MyStudy", "Treatment_vs_Control", "protein")
#' }
#'
#' @export
write_contrast_tables <- function(result_dt, output_dir, project,
                                  contrast, level = "protein") {
  fs::dir_create(output_dir)
  ts  <- Sys.time()

  paths <- c(
    all  = fs::path(output_dir,
                    safe_filename(project,"STAT",level,contrast,"all_results","tsv",ts)),
    sig  = fs::path(output_dir,
                    safe_filename(project,"STAT",level,contrast,"significant","tsv",ts)),
    up   = fs::path(output_dir,
                    safe_filename(project,"STAT",level,contrast,"upregulated","tsv",ts)),
    down = fs::path(output_dir,
                    safe_filename(project,"STAT",level,contrast,"downregulated","tsv",ts))
  )

  data.table::fwrite(result_dt, paths["all"], sep = "\t")
  data.table::fwrite(result_dt[significant == TRUE], paths["sig"], sep = "\t")
  data.table::fwrite(result_dt[direction == "up"],   paths["up"],  sep = "\t")
  data.table::fwrite(result_dt[direction == "down"], paths["down"],sep = "\t")

  logger::log_info("Wrote contrast tables for {contrast} to {output_dir}")
  invisible(paths)
}

#' Build a GSEA-ranked statistic from a result table
#'
#' Computes `sign(log2FC) × -log10(P.Value)` as the ranking metric. Named
#' by `gene_symbol` (or `feature_id` if gene symbol is missing). Ties are
#' broken by fold change magnitude.
#'
#' @param result_dt `data.table` from `run_limma_deqms()`.
#'
#' @return A named numeric vector, sorted descending.
#'
#' @examples
#' dt <- data.table::data.table(
#'   feature_id  = c("P1","P2","P3"),
#'   gene_symbol = c("BRCA1","TP53","EGFR"),
#'   log2FC      = c(2, -1.5, 0.5),
#'   P.Value     = c(0.001, 0.01, 0.5)
#' )
#' build_gsea_ranking(dt)
#'
#' @export
build_gsea_ranking <- function(result_dt) {
  dt <- data.table::copy(result_dt)

  id_col <- if ("gene_symbol" %in% names(dt) &&
                any(!is.na(dt$gene_symbol))) "gene_symbol" else "feature_id"

  dt[, rank_stat := sign(log2FC) * (-log10(pmax(P.Value, 1e-300)))]
  dt <- dt[!is.na(rank_stat)]
  dt <- dt[!duplicated(dt[[id_col]])]
  data.table::setorder(dt, -rank_stat)

  stats::setNames(dt$rank_stat, dt[[id_col]])
}

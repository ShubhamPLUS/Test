#' Export a complete analysis results Excel workbook
#'
#' Writes the canonical `complete_analysis_results.xlsx` with frozen headers,
#' autofilters, and adjusted column widths. Every required sheet is included;
#' missing data results in an empty (but present) sheet rather than omission.
#'
#' @param pfd          A `ProteoForgeData` object.
#' @param stat_results Named list of contrast result data.tables (from
#'   `run_limma_deqms()`). Can be `NULL`.
#' @param qc_results   Named list from `compute_qc_metrics()`. Can be `NULL`.
#' @param enrichment   Named list of enrichment results. Can be `NULL`.
#' @param output_path  Full path for the `.xlsx` file.
#' @param config       ProteoForge config list.
#'
#' @return Invisible path to the written workbook.
#'
#' @examples
#' \dontrun{
#' export_excel(pfd, stat_results, qc_results, NULL,
#'              "results/10_excel/complete_analysis_results.xlsx")
#' }
#'
#' @export
export_excel <- function(pfd, stat_results = NULL, qc_results = NULL,
                         enrichment = NULL, output_path, config = list()) {
  if (!requireNamespace("openxlsx2", quietly = TRUE)) {
    stop("Package 'openxlsx2' required.", call. = FALSE)
  }

  fs::dir_create(fs::path_dir(output_path))
  wb <- openxlsx2::wb_workbook()

  # ── Helper: add styled sheet ──────────────────────────────────────────────
  .add_sheet <- function(wb, sheet_name, data, ...) {
    nm <- substr(sheet_name, 1, 31)  # Excel sheet name limit
    wb$add_worksheet(nm)
    if (!is.null(data) && nrow(data) > 0) {
      wb$add_data_table(nm, x = as.data.frame(data),
                        table_style = "TableStyleLight2",
                        with_filter = TRUE, ...)
    } else {
      wb$add_data(nm, x = data.frame(note = "No data for this sheet"))
    }
    wb
  }

  # ── 1. README ─────────────────────────────────────────────────────────────
  readme <- data.frame(
    Key   = c("Package","Version","Run date","Project","Analysis level",
              "Source software","Parameter hash"),
    Value = c("proteoforge",
              tryCatch(as.character(utils::packageVersion("proteoforge")),
                       error=function(e) "dev"),
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              config$project$name %||% "unknown",
              pfd@analysis_level,
              pfd@source_software,
              pfd@parameter_hash)
  )
  wb <- .add_sheet(wb, "README", readme)

  # ── 2. Sample metadata ─────────────────────────────────────────────────────
  wb <- .add_sheet(wb, "Sample_Metadata", pfd@sample_metadata)

  # ── 3. Feature metadata ────────────────────────────────────────────────────
  wb <- .add_sheet(wb, "Feature_Metadata", pfd@feature_metadata)

  # ── 4. Validation report ───────────────────────────────────────────────────
  val <- pfd@validation
  val_dt <- data.table::data.table(
    field = c("passed","n_samples","n_conditions"),
    value = c(as.character(val$passed %||% ""),
              as.character(val$checks$n_samples %||% ""),
              as.character(val$checks$n_conditions %||% ""))
  )
  wb <- .add_sheet(wb, "Validation_Report", val_dt)

  # ── 5. QC tables ───────────────────────────────────────────────────────────
  if (!is.null(qc_results)) {
    wb <- .add_sheet(wb, "QC_Sample_Summary", qc_results$id_counts)
    wb <- .add_sheet(wb, "Missingness_By_Sample",
                     qc_results$missingness$by_sample)
    wb <- .add_sheet(wb, "Missingness_By_Cond",
                     qc_results$missingness$by_condition)
    wb <- .add_sheet(wb, "QC_Outlier_Report", qc_results$outlier_report)
  }

  # ── 6. Matrices ────────────────────────────────────────────────────────────
  .mat_to_dt <- function(m) {
    dt <- data.table::as.data.table(m, keep.rownames = "feature_id")
    dt
  }
  wb <- .add_sheet(wb, "Raw_Matrix",
                   if (!is.null(pfd@parameters$pre_norm_matrix))
                     .mat_to_dt(pfd@parameters$pre_norm_matrix)
                   else .mat_to_dt(pfd@raw_matrix))
  wb <- .add_sheet(wb, "Normalized_Matrix",
                   if (!is.null(pfd@parameters$pre_impute_matrix))
                     .mat_to_dt(pfd@parameters$pre_impute_matrix)
                   else .mat_to_dt(pfd@raw_matrix))
  wb <- .add_sheet(wb, "Imputed_Matrix", .mat_to_dt(pfd@raw_matrix))

  # Imputation summary
  imp_stats <- pfd@parameters$imputation_stats
  if (!is.null(imp_stats)) {
    imp_dt <- data.table::data.table(
      metric = names(imp_stats),
      value  = as.character(unlist(imp_stats))
    )
    wb <- .add_sheet(wb, "Imputation_Summary", imp_dt)
  }

  # ── 7. Statistical results ─────────────────────────────────────────────────
  if (!is.null(stat_results)) {
    for (cname in names(stat_results)) {
      clean <- sanitise_name(cname)
      wb <- .add_sheet(wb,
                       substr(paste0("DE_All_",  clean), 1, 31),
                       stat_results[[cname]])
      wb <- .add_sheet(wb,
                       substr(paste0("DE_Sig_",  clean), 1, 31),
                       stat_results[[cname]][significant == TRUE])
      wb <- .add_sheet(wb,
                       substr(paste0("DE_Up_",   clean), 1, 31),
                       stat_results[[cname]][direction == "up"])
      wb <- .add_sheet(wb,
                       substr(paste0("DE_Down_", clean), 1, 31),
                       stat_results[[cname]][direction == "down"])
    }
  }

  # ── 8. Enrichment ─────────────────────────────────────────────────────────
  if (!is.null(enrichment)) {
    for (cname in names(enrichment)) {
      clean <- sanitise_name(cname)
      if (!is.null(enrichment[[cname]]$ora)) {
        wb <- .add_sheet(wb, substr(paste0("ORA_",  clean), 1, 31),
                         enrichment[[cname]]$ora)
      }
      if (!is.null(enrichment[[cname]]$gsea)) {
        wb <- .add_sheet(wb, substr(paste0("GSEA_", clean), 1, 31),
                         enrichment[[cname]]$gsea)
      }
    }
  }

  # ── 9. Warnings & errors ───────────────────────────────────────────────────
  warn_dt <- data.table::data.table(message = pf_env()$warnings)
  err_dt  <- data.table::data.table(message = pf_env()$errors)
  wb <- .add_sheet(wb, "Warnings", warn_dt)
  wb <- .add_sheet(wb, "Errors",   err_dt)

  # ── Save ──────────────────────────────────────────────────────────────────
  wb$save(output_path)
  logger::log_info("Excel workbook written to {output_path}")
  invisible(output_path)
}

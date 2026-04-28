#' Run the full ProteoForge analysis pipeline
#'
#' Top-level orchestrator that calls each pipeline module in order.
#' This is the function called by `scripts/run_pipeline.R` and
#' `scripts/run_demo.R`.
#'
#' @param config      Named list from `yaml::read_yaml()` following the
#'   ProteoForge config schema.
#' @param resume      Logical. Resume from a prior targets cache.
#' @param report_only Logical. Only regenerate reports from cached results.
#' @param params_hash Optional sha256 hash to verify reproducibility against.
#'
#' @return Invisible path to the run output directory.
#'
#' @examples
#' \dontrun{
#' cfg <- yaml::read_yaml(
#'   system.file("config","default_lfq.yml", package="proteoforge")
#' )
#' cfg$input$quant_file <- system.file(
#'   "extdata","generic_lfq_small","protein_matrix.tsv",
#'   package="proteoforge")
#' cfg$input$metadata_file <- system.file(
#'   "extdata","generic_lfq_small","sample_metadata.tsv",
#'   package="proteoforge")
#' run_pipeline(cfg)
#' }
#'
#' @export
run_pipeline <- function(config,
                         resume      = FALSE,
                         report_only = FALSE,
                         params_hash = NULL) {

  # ── 0. Initialise ──────────────────────────────────────────────────────────
  seed <- config$project$seed %||% 1234L
  set.seed(seed)

  # Reset accumulated state
  pf_env()$warnings <- character(0)
  pf_env()$errors   <- character(0)

  ts     <- Sys.time()
  dirs   <- make_output_dirs(
    base_dir  = config$project$output_dir %||% "results",
    project   = config$project$name %||% "proteoforge",
    timestamp = ts
  )
  init_logger(dirs["logs"])
  logger::log_info("ProteoForge pipeline starting")
  logger::log_info("Project: {config$project$name}")
  logger::log_info("Output:  {dirs['root']}")

  run_hash <- hash_params(config)
  write_param_hash(run_hash, dirs["root"])
  logger::log_info("Parameter hash: {run_hash}")

  if (!is.null(params_hash) && params_hash != run_hash) {
    stop(sprintf("Parameter hash mismatch. Expected %s, got %s.",
                 params_hash, run_hash), call. = FALSE)
  }

  yaml::write_yaml(config, fs::path(dirs["inputs"], "config_used.yml"))

  # ── 1. Import ──────────────────────────────────────────────────────────────
  logger::log_info("Step 1: Importing data")
  pfd <- .import_dispatch(config)

  # ── 2. Validate metadata ───────────────────────────────────────────────────
  logger::log_info("Step 2: Validating metadata")
  pfd <- validate_sample_metadata(pfd, output_dir = dirs["inputs"])

  stat_results  <- NULL
  qc_results    <- NULL
  all_plots     <- list()
  enrich_results<- NULL

  if (!report_only) {

    # ── 3. Preprocess ────────────────────────────────────────────────────────
    logger::log_info("Step 3: Preprocessing")
    pfd <- filter_features(pfd, config)
    pfd <- log2_transform(pfd)
    pfd <- normalise_matrix(pfd,
                            method = config$preprocessing$normalization %||% "median")
    pfd <- impute_missing(pfd,
                          method = config$preprocessing$imputation %||% "mixed",
                          seed   = seed)

    # Save matrices
    data.table::fwrite(
      data.table::as.data.table(pfd@raw_matrix, keep.rownames="feature_id"),
      fs::path(dirs["prep_tables"], "imputed_matrix.tsv"), sep="\t")

    # ── 4. QC ────────────────────────────────────────────────────────────────
    logger::log_info("Step 4: QC")
    qc_results <- compute_qc_metrics(pfd)

    # Write QC tables
    data.table::fwrite(qc_results$id_counts,
                       fs::path(dirs["qc_tables"],"id_counts.tsv"), sep="\t")
    data.table::fwrite(qc_results$missingness$by_sample,
                       fs::path(dirs["qc_tables"],"missingness_by_sample.tsv"),sep="\t")
    data.table::fwrite(qc_results$missingness$by_condition,
                       fs::path(dirs["qc_tables"],"missingness_by_condition.tsv"),sep="\t")
    data.table::fwrite(qc_results$outlier_report,
                       fs::path(dirs["qc_tables"],"outlier_report.tsv"),sep="\t")

    # QC plots
    all_plots <- generate_qc_plots(pfd, qc_results, dirs["qc_plots"], config)

    # ── 5. Statistics ────────────────────────────────────────────────────────
    if (length(config$statistics$contrasts %||% list()) > 0) {
      logger::log_info("Step 5: Differential abundance")
      stat_results <- run_limma_deqms(pfd, config)

      # Write per-contrast tables and plots
      for (cname in names(stat_results)) {
        write_contrast_tables(stat_results[[cname]],
                              dirs["stat_tables"], config$project$name %||% "proteoforge",
                              cname, pfd@analysis_level)

        vcano <- plot_volcano(stat_results[[cname]],
                              title = sprintf("Volcano: %s", cname),
                              alpha       = config$statistics$alpha      %||% 0.05,
                              lfc_cutoff  = config$statistics$lfc_cutoff %||% 1,
                              label_top_n = config$plots$label_top_n     %||% 15)
        save_plot(vcano, dirs["stat_plots"],
                  safe_filename(config$project$name %||% "proteoforge",
                                "STAT", pfd@analysis_level, cname, "volcano","",ts) |>
                    fs::path_ext_remove(),
                  formats   = config$plots$formats %||% c("pdf","png"),
                  dpi       = config$plots$dpi %||% 300,
                  module    = "STAT", plot_type = "volcano",
                  contrast  = cname, level = pfd@analysis_level,
                  params_hash = run_hash)
        all_plots[[paste0("volcano_", cname)]] <- vcano

        hm <- suppressWarnings(plot_heatmap(pfd, stat_results[[cname]],
                                            top_n = config$plots$heatmap_top_n %||% 50,
                                            title = sprintf("Top DE — %s", cname)))
        if (!is.null(hm)) {
          all_plots[[paste0("heatmap_", cname)]] <- hm
        }

        # GSEA ranking file
        ranked <- build_gsea_ranking(stat_results[[cname]])
        data.table::fwrite(
          data.table::data.table(gene=names(ranked), rank_stat=ranked),
          fs::path(dirs["stat_tables"],
                   safe_filename(config$project$name %||% "proteoforge",
                                 "STAT", pfd@analysis_level, cname, "gsea_ranking","tsv",ts)),
          sep="\t")
      }

      # ── Overlap analysis (multi-contrast) ─────────────────────────────────
      if (length(stat_results) >= 2) {
        for (mode in c("all","up","down")) {
          p_venn <- tryCatch(
            plot_overlap(stat_results, mode=mode,
                         title=sprintf("%s — significant (%s)", config$project$name, mode)),
            error = function(e) NULL
          )
          if (!is.null(p_venn)) {
            save_plot(p_venn, dirs["overlap"],
                      safe_filename(config$project$name %||% "proteoforge",
                                    "OVERLAP", pfd@analysis_level, "multi",
                                    paste0("venn_",mode),"",ts) |> fs::path_ext_remove(),
                      formats=config$plots$formats %||% c("pdf","png"),
                      module="OVERLAP", plot_type=paste0("venn_",mode),
                      params_hash=run_hash)
          }
          int_tbl <- tryCatch(
            get_intersection_table(stat_results, mode=mode),
            error = function(e) NULL
          )
          if (!is.null(int_tbl)) {
            data.table::fwrite(int_tbl,
              fs::path(dirs["overlap"],
                       safe_filename(config$project$name %||% "proteoforge",
                                     "OVERLAP", pfd@analysis_level, "multi",
                                     paste0("intersections_",mode),"tsv",ts)),
              sep="\t")
          }
        }
      }
    }

    # ── 6. Enrichment (stub when packages absent) ─────────────────────────
    if (isTRUE(config$enrichment$enabled)) {
      logger::log_info("Step 6: Functional enrichment")
      enrich_results <- tryCatch(
        run_enrichment(stat_results, pfd, dirs, config),
        error = function(e) {
          collect_warning(sprintf("Enrichment failed: %s", conditionMessage(e)),
                          "enrichment")
          NULL
        }
      )
    }

    # ── 7. Networks ───────────────────────────────────────────────────────
    if (isTRUE(config$network$enabled)) {
      logger::log_info("Step 7: Network analysis")
      tryCatch(
        run_networks(stat_results, pfd, dirs, config),
        error = function(e) collect_warning(sprintf("Network failed: %s",
          conditionMessage(e)), "network")
      )
    }

    # ── 8. Phospho ───────────────────────────────────────────────────────
    if (isTRUE(config$phospho$enabled) ||
        identical(config$input$analysis_level, "phosphosite")) {
      logger::log_info("Step 8: Phosphoproteomics")
      tryCatch(
        run_phospho_pipeline(pfd, stat_results, dirs, config),
        error = function(e) collect_warning(sprintf("Phospho failed: %s",
          conditionMessage(e)), "phospho")
      )
    }
  }

  # ── 9. Reports ──────────────────────────────────────────────────────────────
  logger::log_info("Step 9: Generating reports")

  # Save serialised objects
  saveRDS(pfd, fs::path(dirs["serialized"], "proteoforge_data.rds"))
  if (!is.null(stat_results)) {
    saveRDS(stat_results, fs::path(dirs["serialized"], "stat_results.rds"))
  }

  # Excel
  if (isTRUE(config$reports$excel %||% TRUE)) {
    tryCatch(
      export_excel(pfd, stat_results, qc_results, enrich_results,
                   fs::path(dirs["excel"], "complete_analysis_results.xlsx"),
                   config),
      error = function(e) collect_warning(
        sprintf("Excel export failed: %s", conditionMessage(e)), "export_excel")
    )
  }

  # Word
  if (isTRUE(config$reports$word %||% TRUE)) {
    tryCatch(
      export_word(pfd, stat_results, qc_results, all_plots,
                  fs::path(dirs["reports"], "final_report.docx"), config),
      error = function(e) collect_warning(
        sprintf("Word export failed: %s", conditionMessage(e)), "export_word")
    )
  }

  # PPTX
  if (isTRUE(config$reports$pptx %||% TRUE)) {
    tryCatch(
      export_pptx(pfd, stat_results, qc_results, all_plots,
                  fs::path(dirs["reports"], "final_presentation.pptx"), config),
      error = function(e) collect_warning(
        sprintf("PPTX export failed: %s", conditionMessage(e)), "export_pptx")
    )
  }

  # ── 10. Finalise ────────────────────────────────────────────────────────────
  write_logs(dirs["logs"], pf_env()$warnings, pf_env()$errors)
  write_session_info(dirs["logs"])

  results_summary <- list(
    project       = config$project$name %||% "proteoforge",
    run_time      = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    parameter_hash= run_hash,
    n_features    = nrow(pfd@raw_matrix),
    n_samples     = ncol(pfd@raw_matrix),
    n_contrasts   = length(stat_results %||% list()),
    n_warnings    = length(pf_env()$warnings),
    n_errors      = length(pf_env()$errors),
    output_dir    = dirs["root"]
  )
  jsonlite::write_json(results_summary,
                        fs::path(dirs["root"], "results_summary.json"),
                        pretty=TRUE, auto_unbox=TRUE)

  if (isTRUE(config$reports$zip_results %||% FALSE)) {
    zip_path <- paste0(dirs["root"], ".zip")
    tryCatch(
      utils::zip(zip_path, dirs["root"], flags="-r9q"),
      error = function(e) collect_warning("ZIP failed", "zip")
    )
  }

  logger::log_info("Pipeline complete. Output: {dirs['root']}")
  invisible(dirs["root"])
}

# ── Private dispatch ────────────────────────────────────────────────────────────

#' @keywords internal
.import_dispatch <- function(config) {
  qf  <- config$input$quant_file
  mf  <- config$input$metadata_file
  fmt <- config$input$format %||% "auto"

  if (is.null(qf) || !file.exists(qf)) {
    stop("Quant file not found: ", qf, call.=FALSE)
  }
  if (is.null(mf) || !file.exists(mf)) {
    stop("Metadata file not found: ", mf, call.=FALSE)
  }

  detected <- detect_software(qf, hint = fmt)
  logger::log_info("Detected input format: {detected}")

  switch(detected,
    diann        = import_diann(qf, mf, config),
    spectronaut  = import_spectronaut(qf, mf, config),
    fragpipe     = import_fragpipe(qf, mf, config),
    maxquant     = import_maxquant(qf, mf, config),
    pd           = import_pd(qf, mf, config),
    skyline      = import_skyline(qf, mf, config),
    import_generic(qf, mf, config)
  )
}

# ── Stub functions (implemented in later phases) ───────────────────────────────

# Importers now fully implemented in their own files (Phase 4)

#' @keywords internal
run_enrichment <- function(stat_results, pfd, dirs, config) {
  # Implemented in Phase 5
  NULL
}
#' @keywords internal
run_networks <- function(stat_results, pfd, dirs, config) {
  # Implemented in Phase 5
  NULL
}
#' @keywords internal
run_phospho_pipeline <- function(pfd, stat_results, dirs, config) {
  # Implemented in Phase 6
  NULL
}

# ── Utility ─────────────────────────────────────────────────────────────────────

#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b

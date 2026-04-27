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
#'   system.file("config", "default_lfq.yml", package = "proteoforge")
#' )
#' cfg$input$quant_file <- system.file(
#'   "extdata", "generic_lfq_small", "protein_matrix.tsv",
#'   package = "proteoforge"
#' )
#' cfg$input$metadata_file <- system.file(
#'   "extdata", "generic_lfq_small", "sample_metadata.tsv",
#'   package = "proteoforge"
#' )
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

  # Record parameter hash
  run_hash <- hash_params(config)
  write_param_hash(run_hash, dirs["root"])
  logger::log_info("Parameter hash: {run_hash}")

  if (!is.null(params_hash) && params_hash != run_hash) {
    stop(sprintf(
      "Parameter hash mismatch. Expected %s, got %s.\n",
      params_hash, run_hash
    ), call. = FALSE)
  }

  # Copy config to inputs
  config_out <- fs::path(dirs["inputs"], "config_used.yml")
  yaml::write_yaml(config, config_out)

  # ── 1. Import ──────────────────────────────────────────────────────────────
  logger::log_info("Step 1: Importing data")
  pfd <- import_data(config)

  # ── 2. Validate metadata ───────────────────────────────────────────────────
  logger::log_info("Step 2: Validating metadata")
  pfd <- validate_metadata(pfd)

  if (report_only) {
    logger::log_warn("--report-only: skipping analysis steps")
  } else {

    # ── 3. Preprocess ──────────────────────────────────────────────────────
    logger::log_info("Step 3: Preprocessing")
    pfd <- preprocess_data(pfd, config)

    # ── 4. QC ──────────────────────────────────────────────────────────────
    logger::log_info("Step 4: QC")
    qc_results <- run_qc(pfd, dirs, config)

    # ── 5. Statistics ──────────────────────────────────────────────────────
    logger::log_info("Step 5: Differential abundance")
    stat_results <- run_statistics(pfd, dirs, config)

    # ── 6. Enrichment ──────────────────────────────────────────────────────
    if (isTRUE(config$enrichment$enabled)) {
      logger::log_info("Step 6: Functional enrichment")
      enrich_results <- run_enrichment(stat_results, dirs, config)
    }

    # ── 7. Networks ────────────────────────────────────────────────────────
    if (isTRUE(config$network$enabled)) {
      logger::log_info("Step 7: Network analysis")
      net_results <- run_networks(stat_results, dirs, config)
    }

    # ── 8. Phospho ─────────────────────────────────────────────────────────
    if (isTRUE(config$phospho$enabled) ||
        identical(config$input$analysis_level, "phosphosite")) {
      logger::log_info("Step 8: Phosphoproteomics")
      phospho_results <- run_phospho(pfd, stat_results, dirs, config)
    }

  }

  # ── 9. Reports ─────────────────────────────────────────────────────────────
  logger::log_info("Step 9: Generating reports")
  generate_reports(pfd, dirs, config)

  # ── 10. Finalise ───────────────────────────────────────────────────────────
  write_logs(dirs["logs"],
             warnings = pf_env()$warnings,
             errors   = pf_env()$errors)
  write_session_info(dirs["logs"])

  if (isTRUE(config$reports$zip_results)) {
    zip_path <- fs::path(dirs["root"], "..", paste0(basename(dirs["root"]), ".zip"))
    utils::zip(zip_path, dirs["root"], flags = "-r9q")
    logger::log_info("Results zipped to {zip_path}")
  }

  logger::log_info("Pipeline complete. Output: {dirs['root']}")
  invisible(dirs["root"])
}

# ── Private stub functions (implemented in Phase 2-3) ─────────────────────────

#' @keywords internal
import_data <- function(config) {
  # Stub — implemented in import_generic.R / import_detect.R
  stop("import_data() not yet implemented. (Phase 2)", call. = FALSE)
}

#' @keywords internal
validate_metadata <- function(pfd) {
  stop("validate_metadata() not yet implemented. (Phase 2)", call. = FALSE)
}

#' @keywords internal
preprocess_data <- function(pfd, config) {
  stop("preprocess_data() not yet implemented. (Phase 2)", call. = FALSE)
}

#' @keywords internal
run_qc <- function(pfd, dirs, config) {
  stop("run_qc() not yet implemented. (Phase 2)", call. = FALSE)
}

#' @keywords internal
run_statistics <- function(pfd, dirs, config) {
  stop("run_statistics() not yet implemented. (Phase 3)", call. = FALSE)
}

#' @keywords internal
run_enrichment <- function(stat_results, dirs, config) {
  stop("run_enrichment() not yet implemented. (Phase 5)", call. = FALSE)
}

#' @keywords internal
run_networks <- function(stat_results, dirs, config) {
  stop("run_networks() not yet implemented. (Phase 5)", call. = FALSE)
}

#' @keywords internal
run_phospho <- function(pfd, stat_results, dirs, config) {
  stop("run_phospho() not yet implemented. (Phase 6)", call. = FALSE)
}

#' @keywords internal
generate_reports <- function(pfd, dirs, config) {
  stop("generate_reports() not yet implemented. (Phase 3)", call. = FALSE)
}

# ── Utility ────────────────────────────────────────────────────────────────────

#' Null-coalescing operator
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b

#!/usr/bin/env Rscript
## ProteoForge CLI
##
## Usage:
##   Rscript proteoforge_cli.R --config config.yml [OPTIONS]
##
## Options:
##   --config <path>      Path to YAML config file (required)
##   --resume             Resume a previous run (skip completed steps)
##   --validate-only      Validate inputs and config, then exit
##   --report-only        Regenerate reports from existing results
##   --dry-run            Print planned steps without executing
##   --cores <n>          Number of parallel workers (default: 1)
##   --log-level <level>  Logging verbosity: DEBUG, INFO, WARN, ERROR (default: INFO)
##   --output-dir <path>  Override output directory from config
##   --help               Show this message and exit

suppressPackageStartupMessages({
  library(optparse, quietly = TRUE)
})

# ── Argument parsing ──────────────────────────────────────────────────────
option_list <- list(
  make_option("--config",        type = "character", default = NULL,
              help = "Path to YAML config file [required]"),
  make_option("--resume",        action = "store_true", default = FALSE,
              help = "Resume a previous run (skip completed targets)"),
  make_option("--validate-only", action = "store_true", default = FALSE,
              dest = "validate_only",
              help = "Validate inputs then exit without running analysis"),
  make_option("--report-only",   action = "store_true", default = FALSE,
              dest = "report_only",
              help = "Regenerate reports from pre-existing results"),
  make_option("--dry-run",       action = "store_true", default = FALSE,
              dest = "dry_run",
              help = "Print planned steps without executing"),
  make_option("--cores",         type = "integer",   default = 1L,
              help = "Number of parallel workers [default: 1]"),
  make_option("--log-level",     type = "character", default = "INFO",
              dest = "log_level",
              help = "Log verbosity: DEBUG|INFO|WARN|ERROR [default: INFO]"),
  make_option("--output-dir",    type = "character", default = NULL,
              dest = "output_dir",
              help = "Override output directory from config")
)

parser <- OptionParser(
  usage       = "Rscript %prog --config config.yml [OPTIONS]",
  option_list = option_list,
  description = "ProteoForge: quantitative proteomics analysis platform"
)
opts <- parse_args(parser, positional_arguments = FALSE)

if (is.null(opts$config)) {
  cat("Error: --config is required.\n\n")
  print_help(parser)
  quit(status = 1)
}
if (!file.exists(opts$config)) {
  cat(sprintf("Error: config file not found: %s\n", opts$config))
  quit(status = 1)
}

# ── Load package ─────────────────────────────────────────────────────────
if (!requireNamespace("proteoforge", quietly = TRUE)) {
  # Dev mode: load from source
  pkg_root <- normalizePath(
    file.path(dirname(sys.frame(1)$ofile), "..", ".."), mustWork = FALSE
  )
  if (file.exists(file.path(pkg_root, "DESCRIPTION"))) {
    devtools::load_all(pkg_root, quiet = TRUE)
  } else {
    stop("proteoforge package not installed and source not found.")
  }
} else {
  library(proteoforge, quietly = TRUE)
}

# ── Read config ──────────────────────────────────────────────────────────
config <- yaml::read_yaml(opts$config)

# Apply CLI overrides
if (!is.null(opts$output_dir)) {
  config$output$dir <- opts$output_dir
}

# ── Initialise logger ────────────────────────────────────────────────────
dirs <- make_output_dirs(config)
init_logger(
  log_dir   = dirs$logs_dir,
  log_level = toupper(opts$log_level)
)
logger::log_info("ProteoForge CLI started")
logger::log_info("Config: {opts$config}")
logger::log_info("Flags: resume={opts$resume}  validate_only={opts$validate_only}  dry_run={opts$dry_run}  report_only={opts$report_only}")

# ── --validate-only ───────────────────────────────────────────────────────
if (opts$validate_only) {
  logger::log_info("Running input validation...")
  config_hash <- hash_params(config)
  logger::log_info("Config hash: {config_hash}")

  quant_file    <- config$input$quant_file
  metadata_file <- config$input$metadata_file

  if (!file.exists(quant_file)) {
    logger::log_error("Quant file not found: {quant_file}")
    quit(status = 1)
  }
  if (!file.exists(metadata_file)) {
    logger::log_error("Metadata file not found: {metadata_file}")
    quit(status = 1)
  }

  smd <- data.table::fread(metadata_file)
  tryCatch({
    validate_sample_metadata(smd, config)
    logger::log_info("Sample metadata validation PASSED")
  }, error = function(e) {
    logger::log_error("Sample metadata validation FAILED: {conditionMessage(e)}")
    quit(status = 1)
  })

  logger::log_info("Validation complete. All checks passed.")
  quit(status = 0)
}

# ── --dry-run ─────────────────────────────────────────────────────────────
if (opts$dry_run) {
  cat("DRY RUN — steps that would be executed:\n\n")
  steps <- c(
    "1. Import data",
    "2. Validate sample metadata",
    "3. Compute QC metrics & plots",
    "4. Filter features",
    "5. Log2 transform",
    "6. Normalise (method: %s)" ,
    "7. Classify missingness",
    "8. Impute (method: %s)",
    "9. Build design + contrast matrices",
    "10. Run limma-DEqMS",
    "11. Write contrast tables",
    "12. Generate volcano / heatmap / overlap plots",
    "13. Run ORA / GSEA enrichment",
    "14. Run STRING + OmniPath networks",
    "15. Export Excel / Word / PPTX",
    "16. Render Quarto HTML report",
    "17. Write AI summary JSON",
    "18. Write param hash"
  )
  cat(sprintf(paste(steps, collapse = "\n"),
              config$preprocessing$normalisation %||% "median",
              config$preprocessing$imputation    %||% "mixed"), "\n")
  if (!is.null(config$phospho$enabled) && isTRUE(config$phospho$enabled)) {
    cat("\nPhospho sub-pipeline:\n")
    cat("  P1. Filter Class-1 sites\n")
    cat("  P2. Parent-protein abundance correction\n")
    cat("  P3. Limma-DEqMS on phosphosites\n")
    cat("  P4. Infer kinase activity (KSEA / KEA3 / decoupleR)\n")
    cat("  P5. Phosphomotif enrichment (rmotifx)\n")
    cat("  P6. Annotate sites (PSP + OmniPath)\n")
  }
  quit(status = 0)
}

# ── --report-only ─────────────────────────────────────────────────────────
if (opts$report_only) {
  logger::log_info("--report-only: regenerating reports from existing results")
  results_dir <- dirs$results_dir
  stat_rds    <- file.path(results_dir, "stat_results.rds")
  qc_rds      <- file.path(results_dir, "qc_metrics.rds")

  if (!file.exists(stat_rds) || !file.exists(qc_rds)) {
    logger::log_error("Existing results not found in {results_dir}. Run without --report-only first.")
    quit(status = 1)
  }
  stat_results <- readRDS(stat_rds)
  qc_metrics   <- readRDS(qc_rds)

  render_quarto_report(stat_results, qc_metrics, dirs = dirs, config = config)
  write_ai_summary_json(stat_results, qc_metrics, dirs = dirs, config = config)
  export_word(stat_results, qc_metrics, dirs = dirs, config = config)
  export_pptx(stat_results, dirs = dirs, config = config)
  logger::log_info("Reports regenerated.")
  quit(status = 0)
}

# ── Full run (or resumed run via targets) ─────────────────────────────────
if (opts$resume || requireNamespace("targets", quietly = TRUE)) {
  logger::log_info("Running via targets pipeline (resume={opts$resume})")
  Sys.setenv(PF_CONFIG = opts$config)

  targets_args <- list()
  if (opts$cores > 1L) {
    if (requireNamespace("future", quietly = TRUE)) {
      future::plan(future::multisession, workers = opts$cores)
      targets_args$callr_function <- NULL
    }
  }

  do.call(targets::tar_make, targets_args)

} else {
  logger::log_info("Running via run_pipeline() (targets not installed)")
  run_pipeline(config_file = opts$config)
}

logger::log_info("ProteoForge CLI finished successfully.")
write_logs(dirs)

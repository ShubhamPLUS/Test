#!/usr/bin/env Rscript
## ProteoForge CLI pipeline runner.
## Usage:
##   Rscript scripts/run_pipeline.R --config config.yml [options]

suppressPackageStartupMessages({
  library(proteoforge)
  library(optparse)
})

opt_list <- list(
  make_option("--config",        default = NULL,
              help = "Path to YAML config file (required)"),
  make_option("--quant",         default = NULL,
              help = "Override quant_file in config"),
  make_option("--metadata",      default = NULL,
              help = "Override metadata_file in config"),
  make_option("--output",        default = NULL,
              help = "Override output_dir in config"),
  make_option("--format",        default = NULL,
              help = "Override input format (auto|diann|spectronaut|...)"),
  make_option("--organism",      default = NULL,
              help = "Override organism (human|mouse|rat|...)"),
  make_option("--resume",        action = "store_true", default = FALSE,
              help = "Resume from targets cache"),
  make_option("--validate-only", action = "store_true", default = FALSE,
              dest = "validate_only",
              help = "Only validate inputs, do not run analysis"),
  make_option("--report-only",   action = "store_true", default = FALSE,
              dest = "report_only",
              help = "Only regenerate reports from cached results"),
  make_option("--dry-run",       action = "store_true", default = FALSE,
              dest = "dry_run",
              help = "Print what would be done without running"),
  make_option("--params-hash",   default = NULL,
              dest = "params_hash",
              help = "Verify reproducibility against a prior run hash")
)

opts <- parse_args(OptionParser(option_list = opt_list))

if (is.null(opts$config)) {
  stop("--config is required. Use --help for usage.", call. = FALSE)
}
if (!file.exists(opts$config)) {
  stop("Config file not found: ", opts$config, call. = FALSE)
}

cfg <- yaml::read_yaml(opts$config)

# Apply CLI overrides
if (!is.null(opts$quant))    cfg$input$quant_file      <- opts$quant
if (!is.null(opts$metadata)) cfg$input$metadata_file   <- opts$metadata
if (!is.null(opts$output))   cfg$project$output_dir    <- opts$output
if (!is.null(opts$format))   cfg$input$format          <- opts$format
if (!is.null(opts$organism)) cfg$project$organism      <- opts$organism

if (opts$dry_run) {
  message("DRY RUN: Would run pipeline with config:")
  cat(yaml::as.yaml(cfg))
  quit(save = "no", status = 0)
}

if (opts$validate_only) {
  message("Running validation only...")
  pfd <- import_data(cfg)
  validated <- validate_metadata(pfd)
  message("Validation complete. No errors.")
  quit(save = "no", status = 0)
}

run_pipeline(cfg,
             resume      = opts$resume,
             report_only = opts$report_only,
             params_hash = opts$params_hash)

message("Pipeline complete.")

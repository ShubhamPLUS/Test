#!/usr/bin/env Rscript
## Run the ProteoForge demo pipeline on the built-in generic_lfq_small dataset.
## Usage: Rscript scripts/run_demo.R [--output results]

suppressPackageStartupMessages({
  library(proteoforge)
  library(optparse)
})

opt_list <- list(
  make_option("--output", default = "results",
              help = "Output directory [default: results]")
)
opts <- parse_args(OptionParser(option_list = opt_list))

# Locate built-in demo data
data_dir   <- system.file("extdata", "generic_lfq_small", package = "proteoforge")
quant_file <- file.path(data_dir, "protein_matrix.tsv")
meta_file  <- file.path(data_dir, "sample_metadata.tsv")

if (!file.exists(quant_file)) {
  stop("Demo data not found. Please reinstall the package.", call. = FALSE)
}

message("=== ProteoForge Demo Run ===")
message("Quant file: ", quant_file)
message("Metadata:   ", meta_file)
message("Output:     ", opts$output)

# Load config and override paths
cfg <- yaml::read_yaml(
  system.file("config", "default_lfq.yml", package = "proteoforge")
)
cfg$input$quant_file    <- quant_file
cfg$input$metadata_file <- meta_file
cfg$project$output_dir  <- opts$output
cfg$project$name        <- "ProteoForge_Demo"
cfg$enrichment$enabled  <- FALSE  # skip network/enrichment for demo speed
cfg$network$enabled     <- FALSE

run_pipeline(cfg)

message("\n=== Demo complete. Check ", opts$output, "/ for results. ===")

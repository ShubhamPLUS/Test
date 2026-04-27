#!/usr/bin/env Rscript
## Launch the ProteoForge Shiny application.
## Usage: Rscript scripts/launch_shiny.R [--port 3838] [--host 0.0.0.0]

suppressPackageStartupMessages({
  library(proteoforge)
  library(optparse)
  library(shiny)
})

opt_list <- list(
  make_option("--port", default = 3838L, type = "integer",
              help = "Port to listen on [default: 3838]"),
  make_option("--host", default = "127.0.0.1",
              help = "Host address [default: 127.0.0.1]")
)
opts <- parse_args(OptionParser(option_list = opt_list))

app_dir <- system.file("shiny", package = "proteoforge")
if (!nzchar(app_dir) || !dir.exists(app_dir)) {
  stop("Shiny app directory not found. Please reinstall the package.", call. = FALSE)
}

message(sprintf("Launching ProteoForge Shiny app at http://%s:%d",
                opts$host, opts$port))
shiny::runApp(app_dir, host = opts$host, port = opts$port,
              launch.browser = interactive())

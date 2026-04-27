#' Sanitise a string for use in file names
#'
#' Replaces characters outside `[A-Za-z0-9._-]` with underscores and
#' collapses consecutive underscores.
#'
#' @param x Character vector to sanitise.
#' @return Sanitised character vector.
#'
#' @examples
#' sanitise_name("Treatment vs Control!")
#'
#' @export
sanitise_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

#' Build a structured output file name
#'
#' Format: `<project>__<MODULE>__<level>__<contrast>__<plot_type>__<timestamp>.<ext>`
#'
#' @param project   Project name string.
#' @param module    Analysis module (e.g. `"STAT"`, `"QC"`, `"ENRICH"`).
#' @param level     Analysis level (e.g. `"protein"`, `"phosphosite"`).
#' @param contrast  Contrast name or `"all"`.
#' @param plot_type Short descriptor (e.g. `"volcano"`, `"heatmap"`).
#' @param ext       File extension without dot (e.g. `"pdf"`, `"tsv"`).
#' @param timestamp POSIXct timestamp; defaults to `Sys.time()`.
#'
#' @return A single character string.
#'
#' @examples
#' safe_filename("MyStudy", "STAT", "protein", "Tx_vs_Ctrl", "volcano", "pdf")
#'
#' @export
safe_filename <- function(project, module, level, contrast, plot_type, ext,
                          timestamp = Sys.time()) {
  ts <- format(timestamp, "%Y%m%d_%H%M%S")
  parts <- c(
    sanitise_name(project),
    sanitise_name(module),
    sanitise_name(level),
    sanitise_name(contrast),
    sanitise_name(plot_type),
    ts
  )
  paste0(paste(parts, collapse = "__"), ".", ext)
}

#' Create the standard ProteoForge output directory tree
#'
#' All sub-directories are created with `fs::dir_create()`.
#'
#' @param base_dir Base results directory.
#' @param project  Project name (used as subdirectory prefix).
#' @param timestamp POSIXct timestamp; defaults to `Sys.time()`.
#'
#' @return Named character vector of sub-directory paths (invisibly).
#'
#' @examples
#' \dontrun{
#' dirs <- make_output_dirs(tempdir(), "MyStudy")
#' }
#'
#' @export
make_output_dirs <- function(base_dir, project, timestamp = Sys.time()) {
  ts     <- format(timestamp, "%Y%m%d_%H%M%S")
  run_id <- paste0(sanitise_name(project), "_", ts)
  root   <- fs::path(base_dir, run_id)

  subdirs <- c(
    logs          = "00_logs",
    inputs        = "01_inputs",
    qc_tables     = "02_qc/tables",
    qc_plots      = "02_qc/plots",
    prep_tables   = "03_preprocessing/tables",
    prep_plots    = "03_preprocessing/plots",
    stat_tables   = "04_statistics/tables",
    stat_plots    = "04_statistics/plots",
    overlap       = "05_overlap",
    enrich_tables = "06_enrichment/tables",
    enrich_plots  = "06_enrichment/plots",
    net_tables    = "07_networks/tables",
    net_plots     = "07_networks/plots",
    net_cyto      = "07_networks/cytoscape",
    phospho_tables= "08_phospho/tables",
    phospho_plots = "08_phospho/plots",
    reports       = "09_reports",
    excel         = "10_excel",
    serialized    = "11_serialized"
  )

  paths <- stats::setNames(
    fs::path(root, subdirs),
    names(subdirs)
  )
  fs::dir_create(paths)
  invisible(c(root = root, paths))
}

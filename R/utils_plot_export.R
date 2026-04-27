#' Save a ggplot2 plot in one or more formats and record it in the manifest
#'
#' Every call appends a row to `plot_manifest.tsv` in `output_dir`.
#'
#' @param plot        A ggplot2 (or compatible) plot object.
#' @param output_dir  Directory in which to save plots and manifest.
#' @param filename_base  Base file name without extension (use `safe_filename()`).
#' @param formats     Character vector of extensions: `"pdf"`, `"png"`, `"tiff"`, `"svg"`.
#' @param width       Plot width in inches.
#' @param height      Plot height in inches.
#' @param dpi         Resolution for raster formats.
#' @param module      Module name for manifest (e.g. `"QC"`).
#' @param plot_type   Short plot type label (e.g. `"pca"`).
#' @param contrast    Contrast name or `"all"`.
#' @param level       Analysis level (e.g. `"protein"`).
#' @param caption     Short caption stored in manifest.
#' @param params_hash sha256 parameter hash for the run.
#' @param source_data Optional `data.frame`/`data.table` to save as TSV alongside plot.
#'
#' @return Named character vector of saved file paths (invisibly).
#'
#' @examples
#' \dontrun{
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(hp, mpg)) + ggplot2::geom_point()
#' save_plot(p, tempdir(), "test__QC__protein__all__scatter", c("pdf", "png"))
#' }
#'
#' @export
save_plot <- function(plot,
                      output_dir,
                      filename_base,
                      formats       = c("pdf", "png"),
                      width         = 7,
                      height        = 5,
                      dpi           = 300,
                      module        = "",
                      plot_type     = "",
                      contrast      = "all",
                      level         = "protein",
                      caption       = "",
                      params_hash   = "",
                      source_data   = NULL) {

  fs::dir_create(output_dir)
  saved <- character(0)

  for (fmt in formats) {
    out_path <- fs::path(output_dir, paste0(filename_base, ".", fmt))
    tryCatch({
      if (fmt == "tiff") {
        ggplot2::ggsave(out_path, plot = plot, width = width, height = height,
                        dpi = 600, compression = "lzw", device = "tiff")
      } else {
        ggplot2::ggsave(out_path, plot = plot, width = width, height = height,
                        dpi = dpi, device = fmt)
      }
      saved <- c(saved, stats::setNames(out_path, fmt))
    }, error = function(e) {
      collect_warning(sprintf("Could not save plot %s.%s: %s",
                               filename_base, fmt, conditionMessage(e)),
                      context = "save_plot")
    })
  }

  if (!is.null(source_data)) {
    tsv_path <- fs::path(output_dir, paste0(filename_base, "_data.tsv"))
    data.table::fwrite(data.table::as.data.table(source_data),
                       tsv_path, sep = "\t")
    saved <- c(saved, data = tsv_path)
  }

  .append_plot_manifest(
    output_dir    = output_dir,
    file_paths    = saved[names(saved) != "data"],
    module        = module,
    plot_type     = plot_type,
    contrast      = contrast,
    level         = level,
    width         = width,
    height        = height,
    dpi           = dpi,
    caption       = caption,
    params_hash   = params_hash
  )

  invisible(saved)
}

#' Append rows to plot_manifest.tsv
#'
#' @keywords internal
.append_plot_manifest <- function(output_dir, file_paths, module, plot_type,
                                  contrast, level, width, height, dpi,
                                  caption, params_hash) {
  manifest_path <- fs::path(output_dir, "..", "plot_manifest.tsv")

  rows <- data.table::data.table(
    file_path       = as.character(file_paths),
    module          = module,
    plot_type       = plot_type,
    contrast        = contrast,
    level           = level,
    width           = width,
    height          = height,
    dpi             = dpi,
    format          = tools::file_ext(as.character(file_paths)),
    function_used   = "save_plot",
    caption         = caption,
    parameters_hash = params_hash,
    created_at      = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )

  if (fs::file_exists(manifest_path)) {
    existing <- data.table::fread(manifest_path, sep = "\t")
    combined <- rbind(existing, rows, fill = TRUE)
  } else {
    combined <- rows
  }
  data.table::fwrite(combined, manifest_path, sep = "\t")
  invisible(manifest_path)
}

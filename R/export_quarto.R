#' Render a Quarto HTML report for a ProteoForge analysis
#'
#' Copies the bundled Quarto template into the output directory, writes
#' a `_params.yml` with all analysis artefacts, then calls
#' `quarto::quarto_render()`.  When Quarto is not installed the function
#' returns `NULL` with a warning rather than erroring.
#'
#' @param stat_results       Named list of per-contrast `data.table`s from
#'   [run_limma_deqms()].
#' @param qc_metrics         Output of [compute_qc_metrics()].
#' @param enrichment_results Optional enrichment results list.
#' @param phospho_results    Optional phospho pipeline results list.
#' @param dirs               Named directory list from [make_output_dirs()].
#' @param config             Pipeline config list.
#'
#' @return Path to the rendered HTML file, or `NULL` on failure.
#'
#' @examples
#' \dontrun{
#' render_quarto_report(stat_results, qc_metrics, dirs = dirs, config = cfg)
#' }
#'
#' @export
render_quarto_report <- function(stat_results,
                                  qc_metrics,
                                  enrichment_results = NULL,
                                  phospho_results    = NULL,
                                  dirs,
                                  config             = list()) {
  if (!requireNamespace("quarto", quietly = TRUE)) {
    collect_warning("Package 'quarto' not installed. HTML report skipped.",
                    "render_quarto_report")
    return(NULL)
  }

  template_src <- system.file("templates", "report_template.qmd",
                               package = "proteoforge")
  if (!file.exists(template_src)) {
    collect_warning("Quarto template not found in package. HTML report skipped.",
                    "render_quarto_report")
    return(NULL)
  }

  report_dir  <- dirs$reports_dir %||% dirs$results_dir
  qmd_dest    <- fs::path(report_dir, "proteoforge_report.qmd")
  params_file <- fs::path(report_dir, "_params.yml")

  fs::file_copy(template_src, qmd_dest, overwrite = TRUE)

  # ── Build parameter list ──────────────────────────────────────────────────
  n_sig_total <- sum(vapply(stat_results, function(dt) {
    if ("significant" %in% names(dt)) sum(dt$significant, na.rm = TRUE) else 0L
  }, integer(1)))

  params <- list(
    project_name       = config$project_name %||% "ProteoForge Analysis",
    n_contrasts        = length(stat_results),
    contrast_names     = names(stat_results),
    n_significant      = n_sig_total,
    n_proteins         = if (length(stat_results) > 0) nrow(stat_results[[1]]) else 0L,
    software           = config$input$software %||% "unknown",
    norm_method        = config$preprocessing$normalisation %||% "median",
    imputation_method  = config$preprocessing$imputation %||% "mixed",
    fdr_threshold      = config$stats$fdr %||% 0.05,
    fc_threshold       = config$stats$log2fc_cutoff %||% 1.0,
    has_enrichment     = !is.null(enrichment_results),
    has_phospho        = !is.null(phospho_results),
    analysis_date      = format(Sys.time(), "%Y-%m-%d %H:%M"),
    proteoforge_version= tryCatch(
      as.character(utils::packageVersion("proteoforge")), error = function(e) "dev"
    )
  )

  yaml::write_yaml(list(params = params), params_file)

  out_path <- tryCatch({
    quarto::quarto_render(
      input            = as.character(qmd_dest),
      output_format    = "html",
      execute_params   = params,
      quiet            = TRUE
    )
    html_out <- fs::path_ext_set(qmd_dest, "html")
    if (file.exists(html_out)) {
      logger::log_info("Quarto report rendered: {html_out}")
      as.character(html_out)
    } else NULL
  }, error = function(e) {
    collect_warning(sprintf("Quarto render failed: %s", conditionMessage(e)),
                    "render_quarto_report")
    NULL
  })

  out_path
}

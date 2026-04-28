#' Generate a Word report from analysis results
#'
#' Writes `final_report.docx` via `officer` + `flextable`. Every templated
#' paragraph is keyed in `inst/templates/methods_text.yml` so a future AI
#' module can override them without touching this function.
#'
#' @param pfd          A `ProteoForgeData` object.
#' @param stat_results Named list of contrast result data.tables.
#' @param qc_results   Named list from `compute_qc_metrics()`.
#' @param plot_list    Named list of ggplot2 objects (QC, volcano, heatmap).
#' @param output_path  Full path for the `.docx` file.
#' @param config       ProteoForge config list.
#'
#' @return Invisible path to the written document.
#'
#' @examples
#' \dontrun{
#' export_word(pfd, stat_results, qc_results, plots,
#'             "results/09_reports/final_report.docx")
#' }
#'
#' @export
export_word <- function(pfd, stat_results = NULL, qc_results = NULL,
                        plot_list = NULL, output_path, config = list()) {
  if (!requireNamespace("officer",   quietly = TRUE)) stop("Package 'officer' required.",   call.=FALSE)
  if (!requireNamespace("flextable", quietly = TRUE)) stop("Package 'flextable' required.", call.=FALSE)

  fs::dir_create(fs::path_dir(output_path))

  # Load methods text templates
  tmpl_path <- system.file("templates","methods_text.yml", package="proteoforge")
  tmpl <- if (file.exists(tmpl_path)) yaml::read_yaml(tmpl_path) else list()

  .sub_params <- function(text, params) {
    for (nm in names(params)) {
      text <- gsub(paste0("\\{\\{", nm, "\\}\\}"), params[[nm]], text)
    }
    text
  }

  params <- list(
    project          = config$project$name       %||% "ProteoForge Analysis",
    organism         = config$project$organism   %||% "human",
    normalization    = config$preprocessing$normalization %||% "median",
    imputation       = config$preprocessing$imputation    %||% "mixed",
    engine           = config$statistics$engine  %||% "limma_deqms",
    n_features       = nrow(pfd@raw_matrix),
    n_samples        = ncol(pfd@raw_matrix),
    analysis_date    = format(Sys.Date(), "%B %d, %Y"),
    parameter_hash   = pfd@parameter_hash
  )

  doc <- officer::read_docx()

  # ── Title ─────────────────────────────────────────────────────────────────
  doc <- doc |>
    officer::body_add_par(params$project, style = "heading 1") |>
    officer::body_add_par(
      sprintf("ProteoForge Analysis Report — %s", params$analysis_date),
      style = "Normal"
    ) |>
    officer::body_add_par("", style = "Normal")

  # ── Project metadata ───────────────────────────────────────────────────────
  doc <- officer::body_add_par(doc, "Project Overview", style = "heading 2")
  meta_ft <- flextable::flextable(data.frame(
    Parameter = c("Project","Organism","Analysis level","Source software",
                  "Features","Samples","Parameter hash"),
    Value     = c(params$project, params$organism, pfd@analysis_level,
                  pfd@source_software, params$n_features, params$n_samples,
                  params$parameter_hash)
  ))
  meta_ft <- flextable::autofit(meta_ft)
  doc     <- flextable::body_add_flextable(doc, meta_ft)
  doc     <- officer::body_add_par(doc, "", style="Normal")

  # ── Methods ────────────────────────────────────────────────────────────────
  doc <- officer::body_add_par(doc, "Methods", style="heading 2")
  methods_text <- tmpl$methods %||%
    "## AI-INSERTION-POINT: methods\n\nProtein intensities were processed using ProteoForge (v{{version}}). Log2 transformation was applied, followed by {{normalization}} normalisation. Missing values were imputed using the {{imputation}} strategy. Differential abundance analysis was performed using {{engine}}."
  doc <- officer::body_add_par(doc,
    .sub_params(methods_text, c(params, list(version="0.1.0"))),
    style = "Normal")
  doc <- officer::body_add_par(doc, "", style="Normal")

  # ── Experimental Design ────────────────────────────────────────────────────
  doc <- officer::body_add_par(doc, "Experimental Design", style="heading 2")
  smd_ft <- flextable::flextable(as.data.frame(
    pfd@sample_metadata[, .SD, .SDcols = intersect(
      c("sample_id","condition","batch","replicate","subject_id"),
      names(pfd@sample_metadata))]))
  smd_ft <- flextable::autofit(smd_ft)
  doc    <- flextable::body_add_flextable(doc, smd_ft)
  doc    <- officer::body_add_par(doc, "", style="Normal")

  # ── QC Summary ────────────────────────────────────────────────────────────
  doc <- officer::body_add_par(doc, "Quality Control", style="heading 2")
  qc_text <- tmpl$qc %||%
    "## AI-INSERTION-POINT: qc\n\nQuality control was performed on {{n_samples}} samples. Missingness and intensity distributions were inspected visually."
  doc <- officer::body_add_par(doc, .sub_params(qc_text, params), style="Normal")

  if (!is.null(qc_results)) {
    miss_ft <- flextable::flextable(
      as.data.frame(qc_results$missingness$by_sample[,
        .SD, .SDcols = c("sample_id","condition","n_detected","pct_missing")]))
    miss_ft <- flextable::autofit(miss_ft)
    doc     <- flextable::body_add_flextable(doc, miss_ft)
  }
  doc <- officer::body_add_par(doc, "", style="Normal")

  # ── Differential Abundance ────────────────────────────────────────────────
  if (!is.null(stat_results)) {
    doc <- officer::body_add_par(doc, "Differential Abundance", style="heading 2")
    de_text <- tmpl$differential_abundance %||%
      "## AI-INSERTION-POINT: differential_abundance\n\nDifferential abundance analysis was performed using {{engine}} with Benjamini-Hochberg correction."
    doc <- officer::body_add_par(doc, .sub_params(de_text, params), style="Normal")

    for (cname in names(stat_results)) {
      dt <- stat_results[[cname]]
      n_sig  <- sum(dt$significant, na.rm=TRUE)
      n_up   <- sum(dt$direction == "up",   na.rm=TRUE)
      n_down <- sum(dt$direction == "down", na.rm=TRUE)
      doc    <- officer::body_add_par(doc,
        sprintf("Contrast: %s — %d significant (%d up, %d down)",
                cname, n_sig, n_up, n_down),
        style = "heading 3")

      top_dt <- dt[seq_len(min(20, nrow(dt))),
                   .SD, .SDcols = intersect(
                     c("feature_id","gene_symbol","log2FC",
                       "adj.P.Val","significant","direction"),
                     names(dt))]
      top_ft <- flextable::flextable(as.data.frame(top_dt))
      top_ft <- flextable::autofit(top_ft)
      top_ft <- flextable::colformat_double(top_ft,
                  j = intersect(c("log2FC","adj.P.Val"), names(top_dt)),
                  digits = 3)
      doc    <- flextable::body_add_flextable(doc, top_ft)
      doc    <- officer::body_add_par(doc, "", style="Normal")
    }
  }

  # ── Limitations & Warnings ────────────────────────────────────────────────
  doc <- officer::body_add_par(doc, "Warnings and Limitations", style="heading 2")
  warns <- pf_env()$warnings
  if (length(warns) == 0) warns <- "No warnings recorded."
  for (w in warns) {
    doc <- officer::body_add_par(doc, paste0("• ", w), style="Normal")
  }

  # ── Session Info ──────────────────────────────────────────────────────────
  doc <- officer::body_add_par(doc, "Reproducibility", style="heading 2")
  doc <- officer::body_add_par(doc,
    sprintf("Parameter hash: %s\nR version: %s",
            params$parameter_hash, R.version$version.string),
    style = "Normal")

  print(doc, target = output_path)
  logger::log_info("Word report written to {output_path}")
  invisible(output_path)
}

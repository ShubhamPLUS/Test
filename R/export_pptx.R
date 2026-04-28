#' Generate a PowerPoint presentation from analysis results
#'
#' Writes `final_presentation.pptx` via `officer`. Editable vector graphics
#' are used where `rvg` is available.
#'
#' @param pfd          A `ProteoForgeData` object.
#' @param stat_results Named list of contrast result data.tables.
#' @param qc_results   Named list from `compute_qc_metrics()`.
#' @param plot_list    Named list of ggplot2 objects.
#' @param output_path  Full path for the `.pptx` file.
#' @param config       ProteoForge config list.
#'
#' @return Invisible path to the written file.
#'
#' @examples
#' \dontrun{
#' export_pptx(pfd, stat_results, qc_results, plots,
#'             "results/09_reports/final_presentation.pptx")
#' }
#'
#' @export
export_pptx <- function(pfd, stat_results = NULL, qc_results = NULL,
                        plot_list = NULL, output_path, config = list()) {
  if (!requireNamespace("officer", quietly = TRUE)) {
    stop("Package 'officer' required.", call. = FALSE)
  }

  fs::dir_create(fs::path_dir(output_path))
  has_rvg <- requireNamespace("rvg", quietly = TRUE)

  prs <- officer::read_pptx()

  .add_title_slide <- function(prs, title, subtitle = "") {
    prs <- officer::add_slide(prs, layout = "Title Slide", master = "Office Theme")
    prs <- officer::ph_with(prs, value = title,    location = officer::ph_location_type("ctrTitle"))
    prs <- officer::ph_with(prs, value = subtitle, location = officer::ph_location_type("subTitle"))
    prs
  }

  .add_content_slide <- function(prs, title, body_text) {
    prs <- officer::add_slide(prs, layout = "Title and Content", master = "Office Theme")
    prs <- officer::ph_with(prs, value = title,     location = officer::ph_location_type("title"))
    prs <- officer::ph_with(prs, value = body_text, location = officer::ph_location_type("body"))
    prs
  }

  .add_plot_slide <- function(prs, title, plot_obj,
                               left=0.5, top=1.5, width=9, height=5.5) {
    prs <- officer::add_slide(prs, layout="Title and Content", master="Office Theme")
    prs <- officer::ph_with(prs, value=title,
                             location=officer::ph_location_type("title"))
    if (has_rvg && inherits(plot_obj, "gg")) {
      prs <- officer::ph_with(prs,
               value = rvg::dml(ggobj = plot_obj),
               location = officer::ph_location(left=left, top=top,
                                                width=width, height=height))
    } else {
      tmp <- tempfile(fileext=".png")
      ggplot2::ggsave(tmp, plot_obj, width=width, height=height, dpi=150)
      prs <- officer::ph_with(prs,
               value = officer::external_img(tmp, width=width, height=height),
               location = officer::ph_location(left=left, top=top,
                                                width=width, height=height))
    }
    prs
  }

  project <- config$project$name %||% "ProteoForge Analysis"

  # ── Title slide ───────────────────────────────────────────────────────────
  prs <- .add_title_slide(prs, project,
    sprintf("ProteoForge v0.1.0 | %s | %s",
            pfd@source_software, format(Sys.Date(), "%Y-%m-%d")))

  # ── Pipeline overview ──────────────────────────────────────────────────────
  prs <- .add_content_slide(prs, "Analysis Pipeline",
    paste(c("1. Data import and validation",
            "2. Quality control",
            "3. Preprocessing: filter → log2 → normalize → impute",
            "4. Differential abundance (limma + DEqMS)",
            "5. Functional enrichment (ORA + GSEA)",
            "6. Network analysis (STRING + OmniPath)",
            "7. Report generation"),
          collapse = "\n"))

  # ── Experimental design ────────────────────────────────────────────────────
  cond_summary <- pfd@sample_metadata[, .N, by=condition]
  prs <- .add_content_slide(prs, "Experimental Design",
    paste(apply(cond_summary, 1, function(r) sprintf("%s: n=%s", r[1], r[2])),
          collapse="\n"))

  # ── QC plots ──────────────────────────────────────────────────────────────
  if (!is.null(plot_list$id_counts)) {
    prs <- .add_plot_slide(prs, "QC: Detected Features per Sample",
                           plot_list$id_counts)
  }
  if (!is.null(plot_list$pca)) {
    prs <- .add_plot_slide(prs, "QC: PCA", plot_list$pca)
  }
  if (!is.null(plot_list$missingness)) {
    prs <- .add_plot_slide(prs, "QC: Missingness", plot_list$missingness)
  }

  # ── DE summary slides ──────────────────────────────────────────────────────
  if (!is.null(stat_results)) {
    for (cname in names(stat_results)) {
      dt    <- stat_results[[cname]]
      n_sig <- sum(dt$significant, na.rm=TRUE)
      n_up  <- sum(dt$direction=="up",   na.rm=TRUE)
      n_down<- sum(dt$direction=="down", na.rm=TRUE)
      prs   <- .add_content_slide(prs,
        sprintf("DE: %s", cname),
        sprintf("Significant: %d  |  Up: %d  |  Down: %d\n\nTop hits:\n%s",
                n_sig, n_up, n_down,
                paste(dt[seq_len(min(5,nrow(dt))),
                         ifelse("gene_symbol" %in% names(dt) &
                                  !is.na(gene_symbol), gene_symbol, feature_id)],
                      collapse=", ")))

      if (!is.null(plot_list[[paste0("volcano_", cname)]])) {
        prs <- .add_plot_slide(prs,
          sprintf("Volcano: %s", cname),
          plot_list[[paste0("volcano_", cname)]])
      }
    }
  }

  # ── Warnings ──────────────────────────────────────────────────────────────
  warns <- pf_env()$warnings
  if (length(warns) > 0) {
    prs <- .add_content_slide(prs, "Warnings",
      paste(paste0("• ", warns), collapse="\n"))
  }

  print(prs, target = output_path)
  logger::log_info("PPTX presentation written to {output_path}")
  invisible(output_path)
}

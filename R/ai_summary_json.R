#' Write AI-ready summary JSON
#'
#' Produces a structured JSON file at `results/<run>/ai_summary.json` that
#' aggregates key analysis metrics for downstream AI-assisted interpretation.
#' The schema is intentionally flat and self-describing so that a language
#' model can generate Methods text, Discussion bullets, or figure captions
#' without needing to parse the raw analysis objects.
#'
#' @param stat_results        Named list of per-contrast `data.table`s.
#' @param qc_metrics          Output of [compute_qc_metrics()].
#' @param enrichment_results  Optional enrichment results list.
#' @param phospho_results     Optional phospho pipeline results list.
#' @param config_hash         SHA-256 hash of analysis parameters (character).
#' @param dirs                Named directory list from [make_output_dirs()].
#' @param config              Pipeline config list.
#'
#' @return Path to the written JSON file (invisibly).
#'
#' @examples
#' \dontrun{
#' write_ai_summary_json(stat_results, qc_metrics, dirs = dirs)
#' }
#'
#' @export
write_ai_summary_json <- function(stat_results,
                                   qc_metrics,
                                   enrichment_results = NULL,
                                   phospho_results    = NULL,
                                   config_hash        = NULL,
                                   dirs,
                                   config             = list()) {
  out_dir  <- dirs$results_dir %||% "."
  out_file <- fs::path(out_dir, "ai_summary.json")

  # ── Per-contrast summaries ────────────────────────────────────────────────
  contrast_summaries <- lapply(names(stat_results), function(cname) {
    dt <- stat_results[[cname]]
    sig_dt <- if ("significant" %in% names(dt)) dt[significant == TRUE] else dt[0]
    up_n   <- if ("direction" %in% names(sig_dt)) sum(sig_dt$direction == "up",   na.rm = TRUE) else NA_integer_
    dn_n   <- if ("direction" %in% names(sig_dt)) sum(sig_dt$direction == "down", na.rm = TRUE) else NA_integer_

    top_up <- character(0)
    top_dn <- character(0)
    gene_col <- intersect(c("GeneSymbol","gene_symbol","Gene"), names(dt))[1]
    if (!is.na(gene_col) && "direction" %in% names(dt) && "adj.P.Val" %in% names(dt)) {
      up_sub <- dt[direction == "up"][order(adj.P.Val)][seq_len(min(10, .N))]
      dn_sub <- dt[direction == "down"][order(adj.P.Val)][seq_len(min(10, .N))]
      top_up <- na.omit(up_sub[[gene_col]])
      top_dn <- na.omit(dn_sub[[gene_col]])
    }

    list(
      contrast         = cname,
      n_tested         = nrow(dt),
      n_significant    = nrow(sig_dt),
      n_up             = up_n,
      n_down           = dn_n,
      top_up_proteins  = as.list(top_up),
      top_down_proteins= as.list(top_dn)
    )
  })

  # ── QC summary ───────────────────────────────────────────────────────────
  qc_summary <- list(
    n_samples          = qc_metrics$n_samples   %||% NA_integer_,
    n_features_raw     = qc_metrics$n_features  %||% NA_integer_,
    median_cv_percent  = round(qc_metrics$median_cv * 100, 2) %||% NA_real_,
    outlier_samples    = as.list(qc_metrics$outlier_samples %||% character(0)),
    pca_variance_pc1   = qc_metrics$pca_variance[1] %||% NA_real_,
    pca_variance_pc2   = qc_metrics$pca_variance[2] %||% NA_real_
  )

  # ── Enrichment summary ────────────────────────────────────────────────────
  enrich_summary <- if (!is.null(enrichment_results)) {
    lapply(names(enrichment_results), function(cname) {
      er  <- enrichment_results[[cname]]
      ora <- er$ora
      list(
        contrast       = cname,
        top_go_bp      = .top_terms(ora, "GO_BP", 5L),
        top_go_mf      = .top_terms(ora, "GO_MF", 5L),
        top_kegg       = .top_terms(ora, "KEGG",  5L),
        top_reactome   = .top_terms(ora, "Reactome", 5L)
      )
    })
  } else list()

  # ── Phospho summary ──────────────────────────────────────────────────────
  phospho_summary <- if (!is.null(phospho_results)) {
    list(
      n_sites_tested    = phospho_results$n_sites_tested     %||% NA_integer_,
      n_sites_sig       = phospho_results$n_sites_sig        %||% NA_integer_,
      n_kinases_changed = phospho_results$n_kinases_changed  %||% NA_integer_,
      top_kinases       = as.list(phospho_results$top_kinases %||% character(0))
    )
  } else list()

  # ── AI insertion points ──────────────────────────────────────────────────
  tmpl_file <- system.file("templates", "methods_text.yml", package = "proteoforge")
  ai_prompts <- if (file.exists(tmpl_file)) yaml::read_yaml(tmpl_file) else list()

  # ── Assemble payload ─────────────────────────────────────────────────────
  payload <- list(
    schema_version     = "1.0",
    generated_at       = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    config_hash        = config_hash %||% NA_character_,
    project_name       = config$project_name %||% "ProteoForge Analysis",
    software           = config$input$software %||% "unknown",
    organism           = config$organism %||% "human",
    analysis_params    = list(
      norm_method      = config$preprocessing$normalisation %||% "median",
      imputation_method= config$preprocessing$imputation    %||% "mixed",
      fdr_threshold    = config$stats$fdr                   %||% 0.05,
      fc_threshold     = config$stats$log2fc_cutoff         %||% 1.0
    ),
    qc_summary         = qc_summary,
    contrast_summaries = contrast_summaries,
    enrichment_summary = enrich_summary,
    phospho_summary    = phospho_summary,
    ai_insertion_points= ai_prompts
  )

  jsonlite::write_json(payload, out_file,
                        auto_unbox = TRUE, pretty = TRUE, na = "null")
  logger::log_info("AI summary JSON written: {out_file}")
  invisible(out_file)
}

# ── Internal helpers ───────────────────────────────────────────────────────

#' @keywords internal
.top_terms <- function(ora_list, db, n) {
  if (is.null(ora_list) || !db %in% names(ora_list)) return(list())
  dt <- ora_list[[db]]
  if (is.null(dt) || nrow(dt) == 0) return(list())
  desc_col <- intersect(c("Description","pathway","ID"), names(dt))[1]
  if (is.na(desc_col)) return(list())
  p_col <- intersect(c("p.adjust","p_adjust","pvalue"), names(dt))[1]
  if (!is.na(p_col)) dt <- dt[order(get(p_col))]
  as.list(head(dt[[desc_col]], n))
}

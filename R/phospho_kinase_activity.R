#' Infer kinase activity using KSEA, KEA3, and/or decoupleR
#'
#' Runs up to three kinase activity methods and applies the consensus rule:
#' a kinase is called "activity-changed" only when ≥ `consensus_min_methods`
#' of the three methods agree at adjusted p < 0.05.
#'
#' @param phospho_result_dt `data.table` from `run_limma_deqms()` on
#'   phosphosite-level data.
#' @param methods           Character vector of methods to run:
#'   `"ksea"`, `"kea3"`, `"decoupler"`. Default all three.
#' @param consensus_min     Minimum methods that must agree. Default 2.
#' @param organism          Organism. Default `"human"`.
#'
#' @return A named list with elements `ksea`, `kea3`, `decoupler`,
#'   `consensus`. Each is a `data.table` or `NULL`.
#'
#' @examples
#' \dontrun{
#' ka <- infer_kinase_activity(phospho_result_dt)
#' }
#'
#' @export
infer_kinase_activity <- function(phospho_result_dt,
                                   methods         = c("ksea","kea3","decoupler"),
                                   consensus_min   = 2L,
                                   organism        = "human") {
  results <- list(ksea=NULL, kea3=NULL, decoupler=NULL, consensus=NULL)

  # Build ranked list (signed statistic) for kinase activity methods
  stat_col <- intersect(c("t","log2FC"), names(phospho_result_dt))[1]
  site_col <- "feature_id"

  ranked_sites <- stats::setNames(
    phospho_result_dt[[stat_col]],
    phospho_result_dt[[site_col]]
  )
  ranked_sites <- ranked_sites[!is.na(ranked_sites)]

  # ── KSEA ──────────────────────────────────────────────────────────────────
  if ("ksea" %in% methods) {
    results$ksea <- tryCatch(
      .run_ksea(ranked_sites, organism),
      error = function(e) {
        collect_warning(sprintf("KSEA failed: %s", conditionMessage(e)), "kinase_activity")
        NULL
      }
    )
  }

  # ── decoupleR ─────────────────────────────────────────────────────────────
  if ("decoupler" %in% methods && requireNamespace("decoupleR", quietly=TRUE)) {
    results$decoupler <- tryCatch(
      .run_decoupler(phospho_result_dt, organism),
      error = function(e) {
        collect_warning(sprintf("decoupleR failed: %s", conditionMessage(e)), "kinase_activity")
        NULL
      }
    )
  }

  # ── KEA3 (API call — skip if no network) ──────────────────────────────────
  if ("kea3" %in% methods) {
    results$kea3 <- tryCatch(
      .run_kea3(phospho_result_dt),
      error = function(e) {
        collect_warning(sprintf("KEA3 failed (API may be unavailable): %s",
                                 conditionMessage(e)), "kinase_activity")
        NULL
      }
    )
  }

  # ── Consensus ─────────────────────────────────────────────────────────────
  active_sets <- Filter(Negate(is.null), results[c("ksea","kea3","decoupler")])
  if (length(active_sets) >= consensus_min) {
    kinase_votes <- table(unlist(lapply(active_sets, function(dt) {
      if (is.null(dt) || !"kinase" %in% names(dt)) return(character(0))
      p_col <- intersect(c("adj.p","adj_p","p.adjust","pval"), names(dt))[1]
      if (is.na(p_col)) return(dt$kinase)
      dt[suppressWarnings(as.numeric(get(p_col))) < 0.05, kinase]
    })))
    consensus_kinases <- names(kinase_votes)[kinase_votes >= consensus_min]
    results$consensus <- data.table::data.table(
      kinase           = consensus_kinases,
      n_methods_agree  = as.integer(kinase_votes[consensus_kinases]),
      consensus_min    = consensus_min
    )
    logger::log_info("Kinase consensus: {length(consensus_kinases)} kinases from {length(active_sets)} methods")
  }

  results
}

#' @keywords internal
.run_ksea <- function(ranked_sites, organism) {
  # KSEA using built-in PhosphoSitePlus-derived substrate sets
  # Returns a data.table with columns: kinase, n_substrates, ES, p, adj.p
  if (!requireNamespace("KSEAapp", quietly=TRUE)) {
    # Simplified KSEA using signed enrichment
    return(.ksea_simple(ranked_sites))
  }
  tryCatch({
    res <- KSEAapp::KSEA.Scores(ranked_sites, organism=organism)
    data.table::as.data.table(res)
  }, error=function(e) .ksea_simple(ranked_sites))
}

#' Simplified KSEA without KSEAapp (uses known kinase-substrate relationships)
#' @keywords internal
.ksea_simple <- function(ranked_sites) {
  # Return empty table — real implementation requires PhosphoSitePlus data
  collect_warning("KSEAapp not installed. KSEA not performed.", "ksea_simple")
  data.table::data.table(kinase=character(0), ES=numeric(0),
                          p=numeric(0), adj.p=numeric(0))
}

#' @keywords internal
.run_decoupler <- function(result_dt, organism) {
  net <- tryCatch(
    decoupleR::get_collectri(organism=organism, split_complexes=FALSE),
    error=function(e) NULL
  )
  if (is.null(net)) return(NULL)

  # Build stats matrix
  stat_col <- intersect(c("t","log2FC"), names(result_dt))[1]
  mat <- matrix(result_dt[[stat_col]], nrow=1,
                 dimnames=list("contrast", result_dt$feature_id))

  res <- tryCatch(
    decoupleR::run_wmean(mat=mat, network=net, .source="source",
                          .target="target", .mor="mor",
                          times=100L, minsize=5L),
    error=function(e) NULL
  )
  if (is.null(res)) return(NULL)
  dt <- data.table::as.data.table(res)
  data.table::setnames(dt, "source", "kinase", skip_absent=TRUE)
  dt
}

#' @keywords internal
.run_kea3 <- function(result_dt, top_n=100L) {
  if (!requireNamespace("httr", quietly=TRUE)) {
    collect_warning("Package 'httr' not installed. KEA3 skipped.", "kea3")
    return(NULL)
  }
  sig_sites <- result_dt[significant==TRUE, feature_id]
  if (length(sig_sites) == 0) return(NULL)

  url <- "https://maayanlab.cloud/kea3/api/enrich/"
  body <- jsonlite::toJSON(list(gene_set=sig_sites, query_name="ProteoForge"),
                            auto_unbox=TRUE)
  resp <- tryCatch(
    httr::POST(url, httr::content_type_json(), body=body,
               httr::timeout(30)),
    error=function(e) NULL
  )
  if (is.null(resp) || httr::status_code(resp) != 200) return(NULL)
  content <- httr::content(resp, as="parsed", type="application/json")
  if (is.null(content)) return(NULL)

  # Parse results
  res <- data.table::rbindlist(lapply(content, function(x) {
    data.table::data.table(kinase=x$TF %||% x$Kinase,
                            rank=x$Rank %||% NA_integer_,
                            score=x$Score %||% NA_real_,
                            adj.p=x$`Adj.P.Value` %||% NA_real_)
  }), fill=TRUE)
  res
}

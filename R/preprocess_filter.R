#' Apply feature-level filters to a ProteoForgeData object
#'
#' Filtering order (as specified):
#' 1. Remove decoys / reverse hits.
#' 2. Remove contaminants.
#' 3. Remove "only identified by site" entries.
#' 4. Q-value cutoff.
#' 5. Minimum unique peptide count.
#' 6. Minimum valid fraction in at least one condition.
#'
#' @param pfd    A `ProteoForgeData` object.
#' @param config Named list following the ProteoForge config schema.
#'
#' @return A `ProteoForgeData` object with filtered `@raw_matrix` and
#'   `@feature_metadata`.
#'
#' @examples
#' qf <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
#'                   package="proteoforge")
#' mf <- system.file("extdata","generic_lfq_small","sample_metadata.tsv",
#'                   package="proteoforge")
#' pfd <- import_generic(qf, mf)
#' pfd <- filter_features(pfd)
#'
#' @export
filter_features <- function(pfd, config = list()) {
  stopifnot(methods::is(pfd, "ProteoForgeData"))

  mat  <- pfd@raw_matrix
  fmd  <- pfd@feature_metadata
  smd  <- pfd@sample_metadata
  cfg  <- config$filtering %||% list()

  n_start <- nrow(mat)
  keep    <- rep(TRUE, n_start)

  # ── 1. Decoys / Reverse ───────────────────────────────────────────────────
  if (isTRUE(cfg$remove_reverse_decoys %||% TRUE)) {
    rev_col <- intersect(c("Reverse","Is Decoy","EG.IsDecoy","Decoy"),
                         names(fmd))[1]
    if (!is.na(rev_col)) {
      is_rev <- .is_flagged(fmd[[rev_col]])
      keep   <- keep & !is_rev
      logger::log_info("Filter: removed {sum(is_rev)} reverse/decoy entries")
    }
  }

  # ── 2. Contaminants ───────────────────────────────────────────────────────
  if (isTRUE(cfg$remove_contaminants %||% TRUE)) {
    con_col <- intersect(c("Potential contaminant","Contaminant","Is Contaminant"),
                         names(fmd))[1]
    if (!is.na(con_col)) {
      is_con <- .is_flagged(fmd[[con_col]])
      keep   <- keep & !is_con
      logger::log_info("Filter: removed {sum(is_con)} contaminant entries")
    }
  }

  # ── 3. Only identified by site ────────────────────────────────────────────
  if (isTRUE(cfg$remove_only_identified_by_site %||% TRUE)) {
    ois_col <- intersect(c("Only identified by site"), names(fmd))[1]
    if (!is.na(ois_col)) {
      is_ois <- .is_flagged(fmd[[ois_col]])
      keep   <- keep & !is_ois
      logger::log_info("Filter: removed {sum(is_ois)} only-identified-by-site entries")
    }
  }

  # ── 4. Q-value cutoff ─────────────────────────────────────────────────────
  qval_cutoff <- cfg$qvalue_cutoff %||% 0.01
  qval_col <- intersect(c("q_value","Q.value","Q-value","qvalue"),
                        names(fmd))[1]
  if (!is.na(qval_col)) {
    vals <- suppressWarnings(as.numeric(fmd[[qval_col]]))
    above <- !is.na(vals) & vals > qval_cutoff
    keep  <- keep & !above
    logger::log_info("Filter: removed {sum(above)} entries above q-value {qval_cutoff}")
  }

  # ── 5. Minimum unique peptides ────────────────────────────────────────────
  min_pep <- cfg$min_unique_peptides %||% 1L
  if (min_pep > 0) {
    pep_col <- intersect(c("n_unique_peptides","Unique peptides",
                            "# Unique Peptides","Razor + unique peptides"),
                         names(fmd))[1]
    if (!is.na(pep_col)) {
      n_pep <- suppressWarnings(as.integer(fmd[[pep_col]]))
      below <- !is.na(n_pep) & n_pep < min_pep
      keep  <- keep & !below
      logger::log_info("Filter: removed {sum(below)} entries with < {min_pep} unique peptides")
    }
  }

  # ── 6. Valid values per condition ─────────────────────────────────────────
  min_frac <- cfg$min_valid_fraction_per_condition %||% 0.5
  if (min_frac > 0) {
    valid_filter <- .valid_fraction_filter(mat[keep, , drop=FALSE],
                                           smd, min_frac)
    # keep is currently length n_start; valid_filter is length sum(keep)
    keep_idx     <- which(keep)
    keep[keep_idx[!valid_filter]] <- FALSE
    n_removed_vf <- sum(!valid_filter)
    logger::log_info("Filter: removed {n_removed_vf} entries failing ",
                     ">={round(min_frac*100)}% valid values in all conditions")
  }

  n_kept <- sum(keep)
  logger::log_info("Filter summary: {n_start} → {n_kept} features retained")

  if (n_kept == 0) {
    stop("All features were removed by filtering. Check filter parameters.", call. = FALSE)
  }
  if (n_kept < 10) {
    collect_warning(sprintf("Only %d features remain after filtering. Results may be unreliable.", n_kept),
                    "filter_features")
  }

  pfd@raw_matrix       <- mat[keep, , drop=FALSE]
  pfd@feature_metadata <- fmd[keep, , drop=FALSE]
  if (nrow(pfd@raw_long) > 0) {
    keep_ids <- rownames(pfd@raw_matrix)
    pfd@raw_long <- pfd@raw_long[feature_id %in% keep_ids]
  }
  pfd
}

#' @keywords internal
.is_flagged <- function(x) {
  if (is.logical(x)) return(x)
  if (is.character(x)) return(x == "+" | tolower(x) == "true")
  as.logical(x)
}

#' @keywords internal
.valid_fraction_filter <- function(mat, smd, min_frac) {
  conds <- unique(smd$condition)
  # For each feature, check if at least one condition meets min_frac
  apply(mat, 1, function(row) {
    any(sapply(conds, function(c) {
      s_ids  <- smd[condition == c, sample_id]
      s_ids  <- intersect(s_ids, colnames(mat))
      vals   <- row[s_ids]
      if (length(vals) == 0) return(FALSE)
      sum(!is.na(vals)) / length(vals) >= min_frac
    }))
  })
}

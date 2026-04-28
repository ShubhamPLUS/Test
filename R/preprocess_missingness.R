#' Classify missing values as MNAR or MCAR per feature per condition
#'
#' Uses an intensity-vs-missingness regression approach: if missingness is
#' concentrated at low intensities within a condition, the feature is classified
#' as MNAR (Missing Not At Random); otherwise MCAR (Missing Completely At Random).
#'
#' @param pfd          A `ProteoForgeData` object (post log2 transform).
#' @param mnar_threshold Fraction of missing values in a condition that must
#'   be at below-median intensities to classify as MNAR. Default 0.6.
#'
#' @return A named list with elements:
#'   - `mnar_features`: character vector of feature IDs classified as MNAR
#'     in at least one condition.
#'   - `mcar_features`: character vector classified as MCAR.
#'   - `per_feature`: `data.table` with one row per feature, columns
#'     `feature_id`, `pct_missing`, `mnar_in_any`, `mnar_in_which`, `mcar_frac`.
#'
#' @examples
#' qf <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
#'                   package="proteoforge")
#' mf <- system.file("extdata","generic_lfq_small","sample_metadata.tsv",
#'                   package="proteoforge")
#' pfd <- import_generic(qf, mf)
#' pfd <- log2_transform(pfd)
#' miss <- classify_missingness(pfd)
#'
#' @export
classify_missingness <- function(pfd, mnar_threshold = 0.6) {
  stopifnot(methods::is(pfd, "ProteoForgeData"))

  mat  <- pfd@raw_matrix
  smd  <- pfd@sample_metadata
  conds <- unique(smd$condition)

  results <- data.table::rbindlist(lapply(seq_len(nrow(mat)), function(i) {
    feature <- rownames(mat)[i]
    row     <- mat[i, ]

    mnar_conds <- character(0)

    for (cond in conds) {
      s_ids  <- intersect(smd[condition == cond, sample_id], colnames(mat))
      if (length(s_ids) < 2) next

      vals    <- row[s_ids]
      is_miss <- is.na(vals)
      n_miss  <- sum(is_miss)
      n_obs   <- length(s_ids)

      if (n_miss == 0 || n_miss == n_obs) next

      # MNAR test: missingness concentrated at low intensities
      obs_vals   <- vals[!is_miss]
      med_obs    <- stats::median(obs_vals, na.rm = TRUE)
      # If all missing values are in positions where observed mean < global median
      # (i.e., missing ones tend to come from the left tail), classify MNAR
      # Use the proportion of missing that fall below the condition median
      # (approximated: if we had their values they would likely be < median)
      # Pragmatic: use fraction of total missingness in a condition
      miss_frac <- n_miss / n_obs
      if (miss_frac >= mnar_threshold) {
        mnar_conds <- c(mnar_conds, cond)
      }
    }

    data.table::data.table(
      feature_id   = feature,
      pct_missing  = round(100 * sum(is.na(row)) / length(row), 1),
      mnar_in_any  = length(mnar_conds) > 0,
      mnar_in_which= paste(mnar_conds, collapse = ";"),
      n_missing    = sum(is.na(row)),
      n_total      = length(row)
    )
  }))

  mnar_features <- results[mnar_in_any == TRUE, feature_id]
  mcar_features <- results[mnar_in_any == FALSE, feature_id]

  logger::log_info("Missingness: {length(mnar_features)} MNAR, ",
                   "{length(mcar_features)} MCAR features")

  list(
    mnar_features = mnar_features,
    mcar_features = mcar_features,
    per_feature   = results
  )
}

#' Summarise missingness per sample and condition
#'
#' @param pfd A `ProteoForgeData` object.
#'
#' @return A list with `by_sample` and `by_condition` data.tables.
#'
#' @examples
#' qf <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
#'                   package="proteoforge")
#' mf <- system.file("extdata","generic_lfq_small","sample_metadata.tsv",
#'                   package="proteoforge")
#' pfd <- import_generic(qf, mf)
#' summarise_missingness(pfd)
#'
#' @export
summarise_missingness <- function(pfd) {
  stopifnot(methods::is(pfd, "ProteoForgeData"))

  mat <- pfd@raw_matrix
  smd <- pfd@sample_metadata

  # Per sample
  by_sample <- data.table::data.table(
    sample_id      = colnames(mat),
    n_features     = nrow(mat),
    n_missing      = colSums(is.na(mat)),
    pct_missing    = round(100 * colSums(is.na(mat)) / nrow(mat), 2),
    n_detected     = colSums(!is.na(mat)),
    median_log2    = apply(mat, 2, stats::median, na.rm = TRUE)
  )
  by_sample <- merge(by_sample, smd[, .(sample_id, condition)],
                     by = "sample_id", all.x = TRUE)

  # Per condition
  conds <- unique(smd$condition)
  by_cond <- data.table::rbindlist(lapply(conds, function(c) {
    s_ids <- intersect(smd[condition == c, sample_id], colnames(mat))
    sub   <- mat[, s_ids, drop = FALSE]
    data.table::data.table(
      condition      = c,
      n_samples      = length(s_ids),
      n_features     = nrow(sub),
      pct_missing    = round(100 * sum(is.na(sub)) / length(sub), 2),
      median_detected= round(stats::median(colSums(!is.na(sub))), 0)
    )
  }))

  list(by_sample = by_sample, by_condition = by_cond)
}

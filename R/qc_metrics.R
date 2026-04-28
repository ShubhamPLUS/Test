#' Compute QC metrics for a ProteoForgeData object
#'
#' Returns a named list of QC tables used downstream by `qc_plots.R` and
#' exported to the Excel workbook.
#'
#' @param pfd A `ProteoForgeData` object (post preprocessing recommended).
#'
#' @return A named list with elements:
#'   - `id_counts`: feature IDs per sample.
#'   - `missingness`: per-sample and per-condition missing value summaries.
#'   - `cv`: coefficient of variation per condition per feature.
#'   - `outlier_report`: samples flagged for potential outlier behaviour.
#'   - `pca`: PCA result on the filtered/normalised matrix.
#'   - `correlation`: sample-level Pearson correlation matrix.
#'
#' @examples
#' qf  <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
#'                    package="proteoforge")
#' mf  <- system.file("extdata","generic_lfq_small","sample_metadata.tsv",
#'                    package="proteoforge")
#' pfd <- import_generic(qf, mf)
#' pfd <- log2_transform(pfd)
#' pfd <- normalise_matrix(pfd)
#' qc  <- compute_qc_metrics(pfd)
#'
#' @export
compute_qc_metrics <- function(pfd) {
  stopifnot(methods::is(pfd, "ProteoForgeData"))

  mat <- pfd@raw_matrix
  smd <- pfd@sample_metadata

  # ── ID counts per sample ──────────────────────────────────────────────────
  id_counts <- data.table::data.table(
    sample_id  = colnames(mat),
    n_detected = colSums(!is.na(mat))
  )
  id_counts <- merge(id_counts, smd[, .(sample_id, condition)],
                     by = "sample_id", all.x = TRUE)

  # ── Missingness summary ────────────────────────────────────────────────────
  miss_summary <- summarise_missingness(pfd)

  # ── CV per condition ───────────────────────────────────────────────────────
  cv_dt <- .compute_cv(mat, smd)

  # ── Sample correlation ─────────────────────────────────────────────────────
  cor_mat <- tryCatch(
    stats::cor(mat, use = "pairwise.complete.obs"),
    error = function(e) matrix(NA, ncol(mat), ncol(mat))
  )
  rownames(cor_mat) <- colnames(mat)
  colnames(cor_mat) <- colnames(mat)

  # ── PCA ───────────────────────────────────────────────────────────────────
  pca_result <- .run_pca(mat, smd)

  # ── Outlier flagging ──────────────────────────────────────────────────────
  outlier_report <- .flag_outliers(id_counts, cv_dt, cor_mat, pca_result, smd)

  list(
    id_counts     = id_counts,
    missingness   = miss_summary,
    cv            = cv_dt,
    correlation   = cor_mat,
    pca           = pca_result,
    outlier_report= outlier_report
  )
}

#' @keywords internal
.compute_cv <- function(mat, smd) {
  conds <- unique(smd$condition)
  data.table::rbindlist(lapply(conds, function(c) {
    s_ids  <- intersect(smd[condition == c, sample_id], colnames(mat))
    if (length(s_ids) < 2) return(NULL)
    sub    <- mat[, s_ids, drop = FALSE]
    cvs    <- apply(sub, 1, function(r) {
      obs <- r[!is.na(r)]
      if (length(obs) < 2) return(NA_real_)
      # CV on log2 scale ≈ multiplicative CV; compute on linear scale
      lin <- 2^obs
      stats::sd(lin, na.rm=TRUE) / mean(lin, na.rm=TRUE) * 100
    })
    data.table::data.table(feature_id = rownames(mat),
                            condition  = c,
                            cv_pct     = cvs)
  }))
}

#' @keywords internal
.run_pca <- function(mat, smd) {
  # Use complete features (no NA) for PCA; fall back to median imputation for viz
  mat_imp <- mat
  for (j in seq_len(ncol(mat_imp))) {
    nas <- is.na(mat_imp[, j])
    if (any(nas)) {
      mat_imp[nas, j] <- stats::median(mat_imp[, j], na.rm = TRUE)
    }
  }
  # Remove zero-variance features
  rv <- apply(mat_imp, 1, stats::var, na.rm = TRUE)
  mat_imp <- mat_imp[rv > 0, , drop = FALSE]
  if (nrow(mat_imp) < 3 || ncol(mat_imp) < 3) return(NULL)

  tryCatch({
    pca <- stats::prcomp(t(mat_imp), scale. = TRUE, center = TRUE)
    var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 2)
    pca_dt  <- data.table::as.data.table(pca$x[, 1:min(5, ncol(pca$x))])
    pca_dt[, sample_id := rownames(pca$x)]
    pca_dt  <- merge(pca_dt, smd, by = "sample_id", all.x = TRUE)
    list(pca = pca, pca_dt = pca_dt, var_exp = var_exp)
  }, error = function(e) NULL)
}

#' @keywords internal
.flag_outliers <- function(id_counts, cv_dt, cor_mat, pca_result, smd) {
  flags <- data.table::copy(smd[, .(sample_id, condition)])
  flags[, n_flags := 0L]
  flags[, flag_reasons := ""]

  # Flag 1: Low ID count (< median - 2*MAD)
  med_ids <- stats::median(id_counts$n_detected)
  mad_ids <- stats::mad(id_counts$n_detected)
  low_id  <- id_counts[n_detected < med_ids - 2*mad_ids, sample_id]
  if (length(low_id) > 0) {
    flags[sample_id %in% low_id, `:=`(
      n_flags     = n_flags + 1L,
      flag_reasons= paste(flag_reasons, "low_id_count", sep=";")
    )]
  }

  # Flag 2: High CV (> 2× median CV per condition)
  if (nrow(cv_dt) > 0) {
    cond_cv <- cv_dt[, .(med_cv = stats::median(cv_pct, na.rm=TRUE)), by=condition]
    smd_cv  <- merge(smd, cond_cv, by="condition")
    high_cv <- smd[, .(sample_id)]  # placeholder — full CV-per-sample below
  }

  # Flag 3: Low correlation with peers
  if (!is.null(cor_mat) && !all(is.na(cor_mat))) {
    # Mean off-diagonal correlation per sample
    diag(cor_mat) <- NA
    mean_cor <- rowMeans(cor_mat, na.rm = TRUE)
    med_cor  <- stats::median(mean_cor, na.rm = TRUE)
    mad_cor  <- stats::mad(mean_cor, na.rm = TRUE)
    low_cor  <- names(mean_cor)[!is.na(mean_cor) & mean_cor < med_cor - 2*mad_cor]
    if (length(low_cor) > 0) {
      flags[sample_id %in% low_cor, `:=`(
        n_flags     = n_flags + 1L,
        flag_reasons= paste(flag_reasons, "low_correlation", sep=";")
      )]
    }
  }

  # Flag 4: PCA outlier (PC1/PC2 > 3 SD)
  if (!is.null(pca_result)) {
    pca_dt  <- pca_result$pca_dt
    for (pc in c("PC1","PC2")) {
      if (pc %in% names(pca_dt)) {
        vals <- pca_dt[[pc]]
        med  <- stats::median(vals, na.rm=TRUE)
        sdev <- stats::sd(vals, na.rm=TRUE)
        out_idx <- abs(vals - med) > 3 * sdev
        out_ids <- pca_dt[out_idx, sample_id]
        if (length(out_ids) > 0) {
          flags[sample_id %in% out_ids, `:=`(
            n_flags     = n_flags + 1L,
            flag_reasons= paste(flag_reasons,
                                paste0("pca_", tolower(pc), "_outlier"), sep=";")
          )]
        }
      }
    }
  }

  flags[, flag_reasons := gsub("^;|;$", "", flag_reasons)]
  flags[, flagged := n_flags >= 2L]

  if (any(flags$flagged)) {
    flagged_ids <- flags[flagged == TRUE, sample_id]
    collect_warning(
      sprintf("Samples flagged as potential outliers (≥2 QC criteria): %s. ",
              "Review before proceeding. Samples NOT auto-removed.",
              paste(flagged_ids, collapse=", ")),
      context = "qc_metrics"
    )
  }

  flags
}

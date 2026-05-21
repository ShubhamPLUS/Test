#' Impute missing values in a ProteoForgeData object
#'
#' Implements a mixed MNAR/MCAR strategy (default), MinProb only, kNN only,
#' or no imputation. The pre-imputation matrix is always preserved.
#'
#' **Mixed strategy** (recommended):
#' - MNAR proteins: MinProb — samples from a downshifted normal distribution
#'   (mean = observed_min - 1.8 × sd, sd = 0.3 × observed_sd) to simulate
#'   proteins near the detection limit.
#' - MCAR proteins: kNN (k=10) weighted by Pearson correlation between samples.
#'
#' A loud warning is issued when > 30% of values in any sample are imputed.
#'
#' @param pfd              A `ProteoForgeData` object (post log2 transform).
#' @param method           One of `"mixed"`, `"minprob"`, `"mindet"`,
#'   `"qrilc"`, `"knn"`, `"rf"`, `"none"`. Default `"mixed"`.
#' @param mnar_threshold   Fraction threshold used to classify MNAR.
#'   Passed to `classify_missingness()`. Default 0.6.
#' @param minprob_downshift Downshift in SDs below mean for MinProb. Default 1.8.
#' @param minprob_sdwidth  SD width for MinProb distribution. Default 0.3.
#' @param knn_k            Number of neighbours for kNN. Default 10.
#' @param seed             Random seed for imputation. Default 1234.
#'
#' @return A `ProteoForgeData` object with imputed `@raw_matrix`.
#'   The pre-imputation matrix is stored in `@parameters$pre_impute_matrix`.
#'   Imputation statistics are stored in `@parameters$imputation_stats`.
#'
#' @examples
#' qf <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
#'                   package="proteoforge")
#' mf <- system.file("extdata","generic_lfq_small","sample_metadata.tsv",
#'                   package="proteoforge")
#' pfd <- import_generic(qf, mf)
#' pfd <- log2_transform(pfd)
#' pfd <- normalise_matrix(pfd)
#' pfd <- impute_missing(pfd)
#'
#' @export
impute_missing <- function(pfd,
                           method            = "mixed",
                           mnar_threshold    = 0.6,
                           minprob_downshift = 1.8,
                           minprob_sdwidth   = 0.3,
                           knn_k             = 10L,
                           seed              = 1234L) {

  stopifnot(methods::is(pfd, "ProteoForgeData"))
  method <- match.arg(method, c("mixed","minprob","mindet","qrilc","knn","rf","none"))

  mat <- pfd@raw_matrix
  pfd@parameters$pre_impute_matrix <- mat

  if (method == "none") {
    logger::log_info("impute_missing: no imputation requested")
    return(pfd)
  }

  n_miss_before <- sum(is.na(mat))
  if (n_miss_before == 0) {
    logger::log_info("impute_missing: no missing values found")
    return(pfd)
  }

  set.seed(seed)

  if (method == "mixed") {
    miss_class <- classify_missingness(pfd, mnar_threshold = mnar_threshold)
    mnar_ids   <- miss_class$mnar_features
    mcar_ids   <- miss_class$mcar_features

    # MNAR → MinProb
    if (length(mnar_ids) > 0) {
      mat[mnar_ids, ] <- .impute_minprob(
        mat[mnar_ids, , drop = FALSE],
        downshift = minprob_downshift,
        sdwidth   = minprob_sdwidth,
        seed      = seed
      )
    }
    # MCAR → kNN
    if (length(mcar_ids) > 0) {
      mat[mcar_ids, ] <- .impute_knn(
        mat[mcar_ids, , drop = FALSE],
        k    = knn_k
      )
    }

  } else if (method == "minprob") {
    mat <- .impute_minprob(mat, downshift = minprob_downshift,
                           sdwidth = minprob_sdwidth, seed = seed)
  } else if (method == "mindet") {
    mat <- .impute_mindet(mat)
  } else if (method == "knn") {
    mat <- .impute_knn(mat, k = knn_k)
  } else if (method == "qrilc") {
    mat <- .impute_qrilc(mat, seed = seed)
  } else if (method == "rf") {
    mat <- .impute_rf(mat, seed = seed)
  }

  n_miss_after   <- sum(is.na(mat))
  n_imputed      <- n_miss_before - n_miss_after
  pct_imputed_by_sample <- round(
    100 * colSums(is.na(pfd@raw_matrix) & !is.na(mat)) / nrow(mat), 2
  )

  # Loud warning for high imputation rate
  if (any(pct_imputed_by_sample > 30)) {
    bad_samps <- names(pct_imputed_by_sample)[pct_imputed_by_sample > 30]
    collect_warning(
      sprintf("More than 30%% of values imputed in samples: %s. ",
              "Results should be interpreted with caution.",
              paste(bad_samps, collapse = ", ")),
      context = "impute_missing"
    )
  }

  logger::log_info("impute_missing: method={method}, imputed {n_imputed} values ",
                   "({round(100*n_imputed/length(mat),1)}%)")

  pfd@raw_matrix <- mat
  pfd@parameters$imputation_method <- method
  pfd@parameters$imputation_stats  <- list(
    n_missing_before      = n_miss_before,
    n_missing_after       = n_miss_after,
    n_imputed             = n_imputed,
    pct_imputed_global    = round(100 * n_imputed / length(mat), 2),
    pct_imputed_by_sample = pct_imputed_by_sample
  )
  pfd
}

# ── Private imputation methods ─────────────────────────────────────────────────

#' @keywords internal
.impute_minprob <- function(mat, downshift = 1.8, sdwidth = 0.3, seed = 1234) {
  set.seed(seed)
  for (j in seq_len(ncol(mat))) {
    col_obs  <- mat[, j]
    obs_vals <- col_obs[!is.na(col_obs)]
    if (length(obs_vals) < 3) next
    mu_imp <- mean(obs_vals) - downshift * stats::sd(obs_vals)
    sd_imp <- sdwidth * stats::sd(obs_vals)
    miss   <- is.na(col_obs)
    col_obs[miss] <- stats::rnorm(sum(miss), mean = mu_imp, sd = sd_imp)
    mat[, j] <- col_obs
  }
  mat
}

#' @keywords internal
.impute_mindet <- function(mat, q = 0.01) {
  for (j in seq_len(ncol(mat))) {
    col_obs <- mat[, j]
    min_val <- stats::quantile(col_obs, q, na.rm = TRUE)
    col_obs[is.na(col_obs)] <- min_val
    mat[, j] <- col_obs
  }
  mat
}

#' @keywords internal
.impute_knn <- function(mat, k = 10L) {
  # Simple column-wise kNN using sample-correlation weights
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)

  if (n_cols < 3) {
    # Fallback to median imputation for tiny matrices
    return(.impute_mindet(mat, q = 0.5))
  }

  # Pre-compute pairwise Pearson correlation between samples
  cor_mat <- tryCatch(
    stats::cor(mat, use = "pairwise.complete.obs"),
    error = function(e) matrix(1, n_cols, n_cols)
  )
  cor_mat[is.na(cor_mat)] <- 0
  diag(cor_mat) <- 0

  for (j in seq_len(n_cols)) {
    miss_rows <- which(is.na(mat[, j]))
    if (length(miss_rows) == 0) next

    # Identify k nearest neighbours by correlation
    neighbours <- order(cor_mat[j, ], decreasing = TRUE)
    neighbours <- neighbours[neighbours != j][seq_len(min(k, n_cols - 1))]

    for (r in miss_rows) {
      neigh_vals <- mat[r, neighbours]
      valid_neigh <- neigh_vals[!is.na(neigh_vals)]
      if (length(valid_neigh) == 0) {
        # Fallback: column median
        mat[r, j] <- stats::median(mat[, j], na.rm = TRUE)
      } else {
        wts <- abs(cor_mat[j, neighbours[!is.na(neigh_vals)]])
        wts[wts == 0] <- 1e-10
        mat[r, j] <- stats::weighted.mean(valid_neigh, wts)
      }
    }
  }
  mat
}

#' @keywords internal
.impute_qrilc <- function(mat, seed = 1234) {
  # Quantile regression imputation for left-censored data
  # Simplified version; for production use imputeLCMD::impute.QRILC
  set.seed(seed)
  for (j in seq_len(ncol(mat))) {
    col <- mat[, j]
    obs <- col[!is.na(col)]
    if (length(obs) < 3) next
    n_miss <- sum(is.na(col))
    if (n_miss == 0) next
    # Sample from the lower tail
    q01 <- stats::quantile(obs, 0.01)
    q10 <- stats::quantile(obs, 0.10)
    col[is.na(col)] <- stats::runif(n_miss, min = q01, max = q10)
    mat[, j] <- col
  }
  mat
}

#' @keywords internal
.impute_rf <- function(mat, seed = 1234) {
  if (!requireNamespace("missForest", quietly = TRUE)) {
    stop("Package 'missForest' required for RF imputation. ",
         "Install with: install.packages('missForest')", call. = FALSE)
  }
  set.seed(seed)
  result <- missForest::missForest(t(mat), verbose = FALSE)
  t(result$ximp)
}

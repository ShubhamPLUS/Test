#' Normalise a log2-transformed intensity matrix
#'
#' Applies one of several normalisation methods to the columns of the intensity
#' matrix stored in a `ProteoForgeData` object. Assumes the matrix has already
#' been log2-transformed.
#'
#' @param pfd    A `ProteoForgeData` object (post log2 transform).
#' @param method One of `"none"`, `"median"`, `"quantile"`, `"vsn"`,
#'   `"cyclic_loess"`. Default `"median"`.
#' @param batch_col Optional column name in `sample_metadata` to use for
#'   visualisation-only batch correction (via `limma::removeBatchEffect`).
#'   Does **not** affect statistical modelling — use as covariate in design
#'   matrix for that purpose.
#'
#' @return A `ProteoForgeData` object with `@raw_matrix` containing normalised
#'   values. The pre-normalisation matrix is stored in `parameters$pre_norm_matrix`.
#'
#' @examples
#' qf <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
#'                   package="proteoforge")
#' mf <- system.file("extdata","generic_lfq_small","sample_metadata.tsv",
#'                   package="proteoforge")
#' pfd <- import_generic(qf, mf)
#' pfd <- log2_transform(pfd)
#' pfd <- normalise_matrix(pfd, method = "median")
#'
#' @export
normalise_matrix <- function(pfd, method = "median", batch_col = NULL) {
  stopifnot(methods::is(pfd, "ProteoForgeData"))
  method <- match.arg(method, c("none","median","quantile","vsn","cyclic_loess"))

  mat <- pfd@raw_matrix
  pfd@parameters$pre_norm_matrix <- mat

  if (method == "none") {
    logger::log_info("normalise_matrix: skipping normalisation (method=none)")
    return(pfd)
  }

  mat_norm <- switch(method,
    median       = .norm_median(mat),
    quantile     = .norm_quantile(mat),
    vsn          = .norm_vsn(mat),
    cyclic_loess = .norm_cyclic_loess(mat)
  )

  logger::log_info("normalise_matrix: applied {method} normalisation")
  pfd@raw_matrix              <- mat_norm
  pfd@parameters$normalization <- method
  pfd
}

#' Log2-transform intensity values
#'
#' Applies `log2(x)` element-wise. Values of 0 or negative are set to `NA`
#' (with a warning if any were present). If the data looks already log2
#' (median > 30 suggests linear scale; median ≤ 30 suggests already log),
#' a warning is issued but the transform is still applied.
#'
#' @param pfd A `ProteoForgeData` object.
#'
#' @return A `ProteoForgeData` object with log2-transformed `@raw_matrix`.
#'
#' @examples
#' qf <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
#'                   package="proteoforge")
#' mf <- system.file("extdata","generic_lfq_small","sample_metadata.tsv",
#'                   package="proteoforge")
#' pfd <- import_generic(qf, mf)
#' pfd <- log2_transform(pfd)
#'
#' @export
log2_transform <- function(pfd) {
  stopifnot(methods::is(pfd, "ProteoForgeData"))

  mat <- pfd@raw_matrix
  med <- stats::median(mat, na.rm = TRUE)

  if (med <= 35 && med > 0) {
    collect_warning(
      sprintf("Median intensity %.1f is low; data may already be log-transformed. ",
              "Applying log2 anyway.", med),
      context = "log2_transform"
    )
  }

  n_non_pos <- sum(mat <= 0, na.rm = TRUE)
  if (n_non_pos > 0) {
    mat[mat <= 0] <- NA
    collect_warning(sprintf("Set %d non-positive values to NA before log2 transform",
                            n_non_pos), context = "log2_transform")
  }

  pfd@raw_matrix <- log2(mat)
  pfd@parameters$log2_transformed <- TRUE
  logger::log_info("log2_transform: applied log2 transform")
  pfd
}

# ── Private normalisation functions ────────────────────────────────────────────

#' @keywords internal
.norm_median <- function(mat) {
  col_medians <- apply(mat, 2, stats::median, na.rm = TRUE)
  global_med  <- stats::median(col_medians, na.rm = TRUE)
  sweep(mat, 2, col_medians - global_med, FUN = "-")
}

#' @keywords internal
.norm_quantile <- function(mat) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' required for quantile normalisation.", call. = FALSE)
  }
  limma::normalizeQuantiles(mat)
}

#' @keywords internal
.norm_vsn <- function(mat) {
  if (!requireNamespace("vsn", quietly = TRUE)) {
    stop("Package 'vsn' required for VSN normalisation. ",
         "Install with BiocManager::install('vsn').", call. = FALSE)
  }
  # vsn expects non-log data — we need to back-transform, normalise, re-log
  collect_warning("VSN normalisation applied to log2 data. ",
                  "Consider using pre-log2 values.", "norm_vsn")
  mat2 <- 2^mat
  mat2[is.na(mat2)] <- NA
  fit <- vsn::vsnMatrix(mat2)
  vsn::exprs(fit)
}

#' @keywords internal
.norm_cyclic_loess <- function(mat) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' required for cyclic loess normalisation.", call. = FALSE)
  }
  limma::normalizeCyclicLoess(mat, method = "fast")
}

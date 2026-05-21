#' Run differential abundance analysis with limma + DEqMS
#'
#' Applies `limma::lmFit` → `limma::contrasts.fit` → `limma::eBayes`, then
#' `DEqMS::spectraCounteBayes` (when peptide counts are available) for
#' peptide-count-weighted variance estimation.
#'
#' @param pfd     A `ProteoForgeData` object (post imputation).
#' @param config  Named list following the ProteoForge config schema.
#'
#' @return A named list, one element per contrast. Each element is a
#'   `data.table` with columns: `feature_id`, `gene_symbol`, `protein_id`,
#'   `log2FC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `B`,
#'   `n_peptides`, `pct_imputed`, `significant`, `direction`.
#'
#' @examples
#' \dontrun{
#' qf  <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
#'                    package="proteoforge")
#' mf  <- system.file("extdata","generic_lfq_small","sample_metadata.tsv",
#'                    package="proteoforge")
#' pfd <- import_generic(qf, mf)
#' pfd <- log2_transform(pfd)
#' pfd <- normalise_matrix(pfd)
#' pfd <- impute_missing(pfd)
#' cfg <- list(statistics = list(
#'   design_formula = "~ 0 + condition",
#'   contrasts = list(list(name="Tx_vs_Ctrl",
#'                         numerator="Treatment",
#'                         denominator="Control")),
#'   p_adjust_method = "BH", alpha = 0.05, lfc_cutoff = 1
#' ))
#' res <- run_limma_deqms(pfd, cfg)
#' }
#'
#' @export
run_limma_deqms <- function(pfd, config = list()) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' required.", call. = FALSE)
  }
  stopifnot(methods::is(pfd, "ProteoForgeData"))

  scfg   <- config$statistics %||% list()
  mat    <- pfd@raw_matrix
  smd    <- pfd@sample_metadata
  fmd    <- pfd@feature_metadata

  # Ensure sample order aligns
  smd    <- smd[match(colnames(mat), smd$sample_id)]

  formula_str   <- scfg$design_formula %||% "~ 0 + condition"
  contrast_defs <- scfg$contrasts      %||% list()
  alpha         <- scfg$alpha          %||% 0.05
  lfc_cutoff    <- scfg$lfc_cutoff     %||% 1
  p_adj_method  <- scfg$p_adjust_method %||% "BH"

  # ── Design matrix ─────────────────────────────────────────────────────────
  design <- build_design_matrix(smd, formula_str)
  logger::log_info("run_limma_deqms: design matrix {nrow(design)}×{ncol(design)}")

  # ── Contrast matrix ───────────────────────────────────────────────────────
  if (length(contrast_defs) == 0) {
    stop("No contrasts defined in config$statistics$contrasts.", call. = FALSE)
  }
  cmat <- build_contrast_matrix(contrast_defs, design)

  # ── Peptide counts for DEqMS ──────────────────────────────────────────────
  pep_col <- intersect(c("n_unique_peptides","Unique peptides",
                          "# Unique Peptides","Razor + unique peptides",
                          "Unique Peptide Count"),
                       names(fmd))[1]
  pep_counts <- if (!is.na(pep_col)) {
    setNames(as.integer(fmd[[pep_col]]), fmd$feature_id)
  } else {
    NULL
  }

  # ── % imputed per feature ─────────────────────────────────────────────────
  pre_imp <- pfd@parameters$pre_impute_matrix
  pct_imp <- if (!is.null(pre_imp)) {
    round(100 * rowSums(is.na(pre_imp[rownames(mat), , drop=FALSE])) /
            ncol(mat), 1)
  } else {
    rep(0, nrow(mat))
  }

  # ── limma fit ─────────────────────────────────────────────────────────────
  fit  <- limma::lmFit(mat, design)
  fit2 <- limma::contrasts.fit(fit, cmat)
  fit3 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)

  # ── DEqMS (if peptide counts available) ───────────────────────────────────
  use_deqms <- !is.null(pep_counts) &&
    requireNamespace("DEqMS", quietly = TRUE)

  if (use_deqms) {
    pep_vec <- pep_counts[rownames(mat)]
    pep_vec[is.na(pep_vec)] <- 1L
    fit3$count <- pep_vec
    tryCatch({
      fit3 <- DEqMS::spectraCounteBayes(fit3)
      logger::log_info("run_limma_deqms: DEqMS variance correction applied")
    }, error = function(e) {
      collect_warning(sprintf("DEqMS failed: %s. Using standard eBayes.",
                               conditionMessage(e)), "limma_deqms")
    })
  } else {
    logger::log_info("run_limma_deqms: using standard limma eBayes (no DEqMS)")
  }

  # ── Extract results per contrast ──────────────────────────────────────────
  results <- vector("list", ncol(cmat))
  names(results) <- colnames(cmat)

  for (cname in colnames(cmat)) {
    tt <- if (use_deqms && !is.null(fit3$sca.t)) {
      DEqMS::outputResult(fit3, coef_col = cname)
    } else {
      as.data.frame(limma::topTable(fit3, coef = cname,
                                     number = Inf, sort.by = "none",
                                     adjust.method = p_adj_method))
    }

    dt <- data.table::as.data.table(tt, keep.rownames = "feature_id")
    # Standardise column names
    setnames_safe <- function(dt, old, new) {
      present <- old[old %in% names(dt)]
      new_n   <- new[old %in% names(dt)]
      if (length(present)) data.table::setnames(dt, present, new_n)
    }
    setnames_safe(dt, c("logFC","log2FC","log.FC"),  rep("log2FC",3))
    setnames_safe(dt, c("P.Value","p.value"),         rep("P.Value",2))
    setnames_safe(dt, c("adj.P.Val","p.adj","padj"),  rep("adj.P.Val",3))

    # Add annotations
    ann_cols <- intersect(c("feature_id","protein_id","gene_symbol",
                             "gene_name","description"),
                          names(fmd))
    if (!"feature_id" %in% names(fmd)) {
      fmd_ann <- data.table::copy(fmd)
      fmd_ann[, feature_id := rownames(pfd@raw_matrix)]
    } else {
      fmd_ann <- fmd[, .SD, .SDcols = ann_cols]
    }
    dt <- merge(dt, fmd_ann, by = "feature_id", all.x = TRUE, sort = FALSE)

    # Add pct_imputed
    dt[, pct_imputed := pct_imp[feature_id]]

    # Add n_peptides
    if (!is.null(pep_counts)) {
      dt[, n_peptides := pep_counts[feature_id]]
    } else {
      dt[, n_peptides := NA_integer_]
    }

    # Significance flags
    if (!"log2FC" %in% names(dt)) dt[, log2FC := get(grep("logFC|log2FC",
                                                             names(dt), value=TRUE)[1])]
    dt[, significant := adj.P.Val < alpha & abs(log2FC) >= lfc_cutoff]
    dt[, direction   := ifelse(!significant, "ns",
                               ifelse(log2FC > 0, "up", "down"))]

    # Sort by significance then fold change
    data.table::setorder(dt, adj.P.Val, -abs(log2FC))

    results[[cname]] <- dt
    n_sig <- sum(dt$significant, na.rm = TRUE)
    logger::log_info("Contrast {cname}: {nrow(dt)} features, {n_sig} significant")
  }

  results
}

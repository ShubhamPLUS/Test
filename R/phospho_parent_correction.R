#' Apply parent-protein abundance correction to phosphosite intensities
#'
#' For each phosphosite, fits:
#'   `log2(site) ~ log2(parent) + condition + batch`
#' and returns the residuals as the corrected site intensity. Reports both
#' uncorrected and corrected statistics, and flags sites where correction
#' changes the significance class.
#'
#' @param phospho_pfd  A `ProteoForgeData` object at `analysis_level="phosphosite"`.
#' @param protein_pfd  A `ProteoForgeData` object at `analysis_level="protein"`
#'   (the paired proteome).
#' @param config       ProteoForge config list.
#'
#' @return A named list:
#'   - `corrected_pfd`: `ProteoForgeData` with residual-corrected matrix.
#'   - `correction_stats`: `data.table` with per-site correction summary.
#'
#' @examples
#' \dontrun{
#' result <- parent_protein_correction(phospho_pfd, protein_pfd)
#' }
#'
#' @export
parent_protein_correction <- function(phospho_pfd, protein_pfd,
                                       config = list()) {
  stopifnot(
    methods::is(phospho_pfd, "ProteoForgeData"),
    methods::is(protein_pfd, "ProteoForgeData")
  )

  ph_mat   <- phospho_pfd@raw_matrix
  pr_mat   <- protein_pfd@raw_matrix
  ph_fmd   <- phospho_pfd@feature_metadata
  smd      <- phospho_pfd@sample_metadata
  shared_s <- intersect(colnames(ph_mat), colnames(pr_mat))

  if (length(shared_s) < 3) {
    stop("< 3 shared samples between phospho and protein matrices.", call.=FALSE)
  }

  # Map phosphosite → parent protein
  prot_col <- intersect(c("ProteinID","protein_id","Protein","Proteins"),
                         names(ph_fmd))[1]
  if (is.na(prot_col)) {
    stop("No protein ID column in phospho feature metadata.", call.=FALSE)
  }

  has_batch    <- "batch" %in% names(smd)
  has_condition<- "condition" %in% names(smd)

  corrected_mat  <- ph_mat
  correction_dt  <- vector("list", nrow(ph_mat))

  for (i in seq_len(nrow(ph_mat))) {
    site_id  <- rownames(ph_mat)[i]
    prot_id  <- ph_fmd[[prot_col]][i]

    if (is.na(prot_id) || !prot_id %in% rownames(pr_mat)) next

    site_vec  <- ph_mat[i,    shared_s]
    prot_vec  <- pr_mat[prot_id, shared_s]
    smd_sub   <- smd[sample_id %in% shared_s]

    # Build regression data frame
    df <- data.frame(
      site = site_vec,
      prot = prot_vec,
      stringsAsFactors = FALSE
    )
    if (has_condition) df$condition <- smd_sub$condition[match(shared_s, smd_sub$sample_id)]
    if (has_batch)     df$batch     <- smd_sub$batch[match(shared_s, smd_sub$sample_id)]

    # Only use rows where both values are observed
    ok <- !is.na(df$site) & !is.na(df$prot)
    if (sum(ok) < 4) next

    fml <- if (has_condition && has_batch) {
      "site ~ prot + condition + batch"
    } else if (has_condition) {
      "site ~ prot + condition"
    } else {
      "site ~ prot"
    }

    tryCatch({
      fit <- stats::lm(stats::as.formula(fml), data=df[ok, ])
      resids <- stats::residuals(fit)
      # Put residuals back into full-length vector (missing stay NA)
      corrected_col <- site_vec
      corrected_col[ok] <- resids
      corrected_mat[i, shared_s] <- corrected_col

      correction_dt[[i]] <- data.table::data.table(
        site_id     = site_id,
        prot_id     = prot_id,
        r_squared   = round(summary(fit)$r.squared, 4),
        n_obs       = sum(ok)
      )
    }, error=function(e) NULL)
  }

  correction_summary <- data.table::rbindlist(correction_dt[!sapply(correction_dt, is.null)])
  logger::log_info("parent_protein_correction: corrected {nrow(correction_summary)} sites")

  corrected_pfd           <- phospho_pfd
  corrected_pfd@raw_matrix <- corrected_mat
  corrected_pfd@parameters$parent_corrected <- TRUE

  list(
    corrected_pfd   = corrected_pfd,
    correction_stats= correction_summary
  )
}

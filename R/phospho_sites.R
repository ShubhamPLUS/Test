#' Filter phosphosites to Class-1 localisation probability
#'
#' Retains only sites with localisation probability ≥ `loc_cutoff` and
#' standardises the site ID to `<UniProt>_<residue><position>` format.
#'
#' @param pfd       A `ProteoForgeData` object with `analysis_level="phosphosite"`.
#' @param loc_cutoff Minimum localisation probability. Default 0.75.
#' @param loc_col    Name of the localisation probability column in
#'   `feature_metadata`. Default auto-detected.
#'
#' @return A `ProteoForgeData` object filtered to Class-1 sites with
#'   standardised site IDs.
#'
#' @examples
#' qf <- system.file("extdata","phospho_lfq_small","phospho_sites.tsv",
#'                   package="proteoforge")
#' mf <- system.file("extdata","phospho_lfq_small","sample_metadata.tsv",
#'                   package="proteoforge")
#' pfd <- import_generic(qf, mf,
#'                       config=list(input=list(analysis_level="phosphosite")),
#'                       feature_id_col="SiteID")
#' pfd <- filter_phospho_sites(pfd, loc_cutoff=0.75)
#'
#' @export
filter_phospho_sites <- function(pfd, loc_cutoff = 0.75, loc_col = NULL) {
  stopifnot(methods::is(pfd, "ProteoForgeData"))

  fmd <- pfd@feature_metadata

  # Auto-detect localisation probability column
  if (is.null(loc_col)) {
    loc_col <- intersect(
      c("LocalizationProb","localization_prob","Localization prob",
        "PTM.Site.Confidence","PTMSiteConfidence"),
      names(fmd)
    )[1]
  }

  if (is.na(loc_col) || !loc_col %in% names(fmd)) {
    collect_warning("No localisation probability column found. Class-1 filter skipped.",
                    "filter_phospho_sites")
    return(pfd)
  }

  loc_vals <- suppressWarnings(as.numeric(fmd[[loc_col]]))
  keep     <- !is.na(loc_vals) & loc_vals >= loc_cutoff
  n_before <- nrow(pfd@raw_matrix)
  n_after  <- sum(keep)

  logger::log_info("filter_phospho_sites: {n_before} → {n_after} Class-1 sites (≥{loc_cutoff})")

  pfd@raw_matrix       <- pfd@raw_matrix[keep, , drop=FALSE]
  pfd@feature_metadata <- fmd[keep]
  if (nrow(pfd@raw_long) > 0) {
    keep_ids <- rownames(pfd@raw_matrix)
    pfd@raw_long <- pfd@raw_long[feature_id %in% keep_ids]
  }
  pfd
}

#' Standardise phosphosite IDs to `<UniProt>_<residue><position>` format
#'
#' @param pfd A `ProteoForgeData` object.
#' @param protein_col  Column in `feature_metadata` with protein IDs.
#' @param residue_col  Column with amino-acid residue (S/T/Y).
#' @param position_col Column with site position integer.
#'
#' @return A `ProteoForgeData` object with standardised `feature_id` values.
#'
#' @examples
#' \dontrun{
#' pfd <- standardise_site_ids(pfd, "ProteinID", "Residue", "Position")
#' }
#'
#' @export
standardise_site_ids <- function(pfd,
                                  protein_col  = "ProteinID",
                                  residue_col  = "Residue",
                                  position_col = "Position") {
  fmd <- pfd@feature_metadata
  req <- c(protein_col, residue_col, position_col)
  missing_c <- setdiff(req, names(fmd))
  if (length(missing_c) > 0) {
    collect_warning(sprintf("standardise_site_ids: columns not found: %s",
                             paste(missing_c, collapse=", ")),
                    "standardise_site_ids")
    return(pfd)
  }

  new_ids <- paste0(fmd[[protein_col]], "_",
                     fmd[[residue_col]],
                     fmd[[position_col]])
  fmd[, feature_id := new_ids]
  rownames(pfd@raw_matrix) <- new_ids

  if ("feature_id" %in% names(pfd@raw_long)) {
    old_ids <- rownames(pfd@raw_matrix)
    pfd@raw_long[, feature_id := new_ids[match(feature_id, old_ids)]]
  }

  pfd@feature_metadata <- fmd
  pfd
}

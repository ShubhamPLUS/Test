#' Overlay OmniPath signaling, complex and kinase-substrate edges
#'
#' Retrieves directed signaling interactions from OmniPath for the set of
#' detected proteins, filters to significant hits, and returns edge/node tables.
#'
#' @param result_dt   `data.table` from `run_limma_deqms()`.
#' @param organism    Organism string. Default `"human"`.
#' @param output_dir  Directory to write files. If `NULL`, no files written.
#'
#' @return A named list: `interactions`, `complexes`, `kinase_substrate`.
#'
#' @examples
#' \dontrun{
#' op <- run_omnipath(result_dt, organism="human")
#' }
#'
#' @export
run_omnipath <- function(result_dt, organism="human", output_dir=NULL) {
  if (!requireNamespace("OmnipathR", quietly=TRUE)) {
    collect_warning("Package 'OmnipathR' not installed. OmniPath skipped.",
                    "run_omnipath")
    return(NULL)
  }

  id_col <- if ("gene_symbol" %in% names(result_dt) &&
                any(!is.na(result_dt$gene_symbol))) "gene_symbol" else "feature_id"
  sig_ids <- unique(na.omit(result_dt[significant==TRUE, get(id_col)]))

  if (length(sig_ids) == 0) return(NULL)

  # ── Signaling interactions ─────────────────────────────────────────────────
  interactions <- tryCatch({
    int <- OmnipathR::import_intercell_network(
      transmitter_categories = c("ligand","receptor"),
      receiver_categories    = c("ligand","receptor")
    )
    int_dt <- data.table::as.data.table(int)
    int_dt <- int_dt[source_genesymbol %in% sig_ids |
                       target_genesymbol %in% sig_ids]
    int_dt
  }, error=function(e) NULL)

  # ── Kinase-substrate ───────────────────────────────────────────────────────
  ks <- tryCatch({
    ks_raw <- OmnipathR::get_kinase_substrate()
    ks_dt  <- data.table::as.data.table(ks_raw)
    ks_dt[enzyme_genesymbol %in% sig_ids | substrate_genesymbol %in% sig_ids]
  }, error=function(e) NULL)

  if (!is.null(output_dir)) {
    fs::dir_create(output_dir)
    if (!is.null(interactions) && nrow(interactions) > 0)
      data.table::fwrite(interactions,
                         fs::path(output_dir,"omnipath_interactions.tsv"), sep="\t")
    if (!is.null(ks) && nrow(ks) > 0)
      data.table::fwrite(ks,
                         fs::path(output_dir,"omnipath_kinase_substrate.tsv"), sep="\t")
  }

  list(
    interactions    = interactions,
    kinase_substrate= ks
  )
}

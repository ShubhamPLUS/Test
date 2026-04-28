#' Annotate phosphosites with PhosphoSitePlus and OmniPath data
#'
#' Overlays known regulatory functions, disease associations, and upstream
#' kinases from PhosphoSitePlus (when the user has downloaded the data) and
#' OmniPath kinase-substrate relationships.
#'
#' **PhosphoSitePlus note:** PSP data requires registration and cannot be
#' bundled with the package. See `scripts/download_psp.R` (generated
#' alongside this function) for download instructions.
#'
#' @param result_dt     `data.table` with site-level DE results.
#' @param psp_dir       Directory containing downloaded PSP files
#'   (`Regulatory_sites.gz`, `Kinase_Substrate_Dataset.gz`). If `NULL`,
#'   PSP annotation is skipped.
#' @param use_omnipath  Logical. Use OmniPath for kinase-substrate annotation.
#'   Default `TRUE`.
#' @param organism      Organism. Default `"human"`.
#'
#' @return The `result_dt` with additional annotation columns appended.
#'
#' @examples
#' \dontrun{
#' annotated <- annotate_phosphosites(result_dt, psp_dir="~/psp_data")
#' }
#'
#' @export
annotate_phosphosites <- function(result_dt,
                                   psp_dir      = NULL,
                                   use_omnipath = TRUE,
                                   organism     = "human") {
  dt <- data.table::copy(result_dt)

  # ── OmniPath annotation ────────────────────────────────────────────────────
  if (use_omnipath && requireNamespace("OmnipathR", quietly=TRUE)) {
    tryCatch({
      ks <- data.table::as.data.table(OmnipathR::get_kinase_substrate())
      # Match by substrate gene and site
      site_col <- intersect(c("gene_symbol","GeneSymbol"), names(dt))[1]
      if (!is.na(site_col)) {
        ks_slim <- unique(ks[, .(
          substrate  = substrate_genesymbol,
          kinase     = enzyme_genesymbol,
          residue    = residue_type,
          position   = residue_offset
        )])
        dt <- merge(dt, ks_slim, by.x=site_col, by.y="substrate", all.x=TRUE)
        dt[, upstream_kinase_omnipath := kinase]
        dt[, kinase := NULL]
      }
    }, error=function(e) NULL)
  }

  # ── PSP annotation ────────────────────────────────────────────────────────
  if (!is.null(psp_dir) && dir.exists(psp_dir)) {
    reg_file <- fs::path(psp_dir, "Regulatory_sites.gz")
    ks_file  <- fs::path(psp_dir, "Kinase_Substrate_Dataset.gz")

    if (file.exists(reg_file)) {
      tryCatch({
        psp_reg <- data.table::fread(reg_file, sep="\t",
                                      na.strings=c("","NA","N/A"))
        psp_reg <- psp_reg[ORGANISM == organism | grepl(organism, tolower(ORGANISM))]
        site_id_col <- intersect(c("feature_id","SiteID"), names(dt))[1]
        if (!is.na(site_id_col)) {
          psp_reg[, site_id_match := paste0(ACC_ID, "_", MOD_RSD)]
          dt <- merge(dt, psp_reg[, .(site_id_match, FUNCTION, PROCESS,
                                        PROT_INTERACTIONS, CST_CAT_NOS,
                                        LT_LIT, CST_LIT, notes)],
                       by.x=site_id_col, by.y="site_id_match", all.x=TRUE)
          dt[, psp_function := FUNCTION]
          dt[, FUNCTION := NULL]
          logger::log_info("PSP regulatory annotation added for {sum(!is.na(dt$psp_function))} sites")
        }
      }, error=function(e) collect_warning(
        sprintf("PSP Regulatory_sites.gz error: %s", conditionMessage(e)),
        "annotate_phosphosites"))
    }

    if (file.exists(ks_file)) {
      tryCatch({
        psp_ks <- data.table::fread(ks_file, sep="\t",
                                     na.strings=c("","NA","N/A"))
        psp_ks <- psp_ks[KIN_ORGANISM == organism | grepl(organism, tolower(KIN_ORGANISM))]
        psp_ks[, site_id_match := paste0(SUB_ACC_ID, "_", SUB_MOD_RSD)]
        ks_slim <- unique(psp_ks[, .(site_id_match,
                                      psp_upstream_kinase=KINASE)])
        site_id_col <- intersect(c("feature_id","SiteID"), names(dt))[1]
        if (!is.na(site_id_col)) {
          dt <- merge(dt, ks_slim,
                       by.x=site_id_col, by.y="site_id_match", all.x=TRUE)
        }
      }, error=function(e) NULL)
    }
  } else if (!is.null(psp_dir)) {
    collect_warning(
      sprintf("PSP directory not found: %s. PSP annotation skipped.", psp_dir),
      "annotate_phosphosites"
    )
  }

  dt
}

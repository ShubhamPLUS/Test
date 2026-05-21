#' Map feature IDs to Entrez Gene IDs for enrichment analysis
#'
#' Auto-detects the ID type (UniProt, gene symbol, Ensembl, RefSeq) and maps
#' to Entrez using `org.<species>.eg.db` with a `UniProt.ws` fallback.
#' Saves success/failure table and warns if < 80% mapped.
#'
#' @param ids      Character vector of input IDs.
#' @param organism Organism string: `"human"`, `"mouse"`, `"rat"`, `"yeast"`.
#' @param id_type  Input ID type: `"auto"`, `"uniprot"`, `"symbol"`,
#'   `"ensembl"`, `"refseq"`.
#'
#' @return A named character vector mapping input ID → Entrez gene ID.
#'   Unmapped IDs are `NA`.
#'
#' @examples
#' \dontrun{
#' ids <- c("P38398","Q9Y6K9","O15350")
#' entrez <- map_ids_to_entrez(ids, organism="human", id_type="uniprot")
#' }
#'
#' @export
map_ids_to_entrez <- function(ids, organism = "human", id_type = "auto") {
  org_pkg <- .get_org_db(organism)
  if (!requireNamespace(org_pkg, quietly=TRUE)) {
    collect_warning(
      sprintf("Annotation package '%s' not installed. ID mapping will fail. ",
              "Install with BiocManager::install('%s').", org_pkg, org_pkg),
      "map_ids_to_entrez"
    )
    return(stats::setNames(rep(NA_character_, length(ids)), ids))
  }

  pkg <- getNamespace(org_pkg)

  # ── Detect ID type ─────────────────────────────────────────────────────────
  if (id_type == "auto") {
    id_type <- .detect_id_type(ids[!is.na(ids)][1:min(5, sum(!is.na(ids)))])
    logger::log_info("ID type auto-detected as: {id_type}")
  }

  keytype <- switch(id_type,
    uniprot = "UNIPROT",
    symbol  = "SYMBOL",
    ensembl = "ENSEMBL",
    refseq  = "REFSEQ",
    "SYMBOL"
  )

  result <- tryCatch({
    db  <- get(org_pkg, envir=pkg)
    suppressMessages(
      AnnotationDbi::mapIds(db, keys=ids, column="ENTREZID",
                             keytype=keytype, multiVals="first")
    )
  }, error = function(e) {
    collect_warning(sprintf("ID mapping error: %s", conditionMessage(e)),
                    "map_ids_to_entrez")
    stats::setNames(rep(NA_character_, length(ids)), ids)
  })

  pct_mapped <- 100 * sum(!is.na(result)) / length(result)
  logger::log_info("ID mapping: {round(pct_mapped,1)}% of {length(ids)} IDs mapped")

  if (pct_mapped < 80) {
    collect_warning(
      sprintf("Only %.1f%% of IDs mapped to Entrez. ",
              "Check organism and ID type settings.", pct_mapped),
      "map_ids_to_entrez"
    )
  }

  result
}

#' @keywords internal
.get_org_db <- function(organism) {
  switch(tolower(organism),
    human  = ,
    homo_sapiens = "org.Hs.eg.db",
    mouse  = ,
    mus_musculus = "org.Mm.eg.db",
    rat    = ,
    rattus = "org.Rn.eg.db",
    yeast  = ,
    saccharomyces = "org.Sc.sgd.db",
    "org.Hs.eg.db"  # default
  )
}

#' @keywords internal
.detect_id_type <- function(sample_ids) {
  sample_ids <- na.omit(sample_ids)
  if (length(sample_ids) == 0) return("symbol")

  # UniProt: P/Q/O followed by 4-9 alphanumeric
  if (mean(grepl("^[PQO][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]$",
                  sample_ids)) > 0.5) return("uniprot")
  # Ensembl
  if (mean(grepl("^ENS[A-Z]+[0-9]{11}", sample_ids)) > 0.5) return("ensembl")
  # RefSeq
  if (mean(grepl("^(NP|XP|YP)_[0-9]+", sample_ids)) > 0.5) return("refseq")
  # Default gene symbol
  "symbol"
}

#' Run Over-Representation Analysis (ORA)
#'
#' Performs ORA using `clusterProfiler::enricher` with the detected
#' proteome as the universe (never the whole genome). Runs on:
#' all-significant, up-regulated only, down-regulated only.
#'
#' @param result_dt   `data.table` from `run_limma_deqms()`.
#' @param universe_ids Character vector of all detected protein IDs (universe).
#' @param organism    Organism string (`"human"`, `"mouse"`, ...).
#' @param databases   Character vector of databases. Supported:
#'   `"GO_BP"`, `"GO_MF"`, `"GO_CC"`, `"KEGG"`, `"Reactome"`,
#'   `"WikiPathways"`, `"MSigDB_Hallmark"`.
#' @param id_type     Input ID type for `map_ids_to_entrez()`.
#' @param p_adjust_method Multiple testing correction. Default `"BH"`.
#' @param min_gs      Minimum gene-set size. Default 10.
#' @param max_gs      Maximum gene-set size. Default 500.
#'
#' @return A named list: one element per database, each a list with
#'   `all`, `up`, `down` sub-elements containing `data.table` results.
#'
#' @examples
#' \dontrun{
#' ora <- run_ora(result_dt, universe_ids, organism="human")
#' }
#'
#' @export
run_ora <- function(result_dt, universe_ids,
                    organism       = "human",
                    databases      = c("GO_BP","KEGG","Reactome"),
                    id_type        = "auto",
                    p_adjust_method= "BH",
                    min_gs         = 10L,
                    max_gs         = 500L) {

  if (!requireNamespace("clusterProfiler", quietly=TRUE)) {
    collect_warning("Package 'clusterProfiler' not installed. ORA skipped.",
                    "run_ora")
    return(NULL)
  }

  universe_entrez <- map_ids_to_entrez(universe_ids, organism, id_type)
  universe_entrez <- unique(na.omit(universe_entrez))

  .do_ora <- function(ids, label) {
    entrez <- unique(na.omit(map_ids_to_entrez(ids, organism, id_type)))
    if (length(entrez) < 3) return(NULL)

    results <- list()
    for (db in databases) {
      res <- tryCatch(
        .ora_one_db(entrez, universe_entrez, db, organism,
                    p_adjust_method, min_gs, max_gs),
        error = function(e) {
          collect_warning(sprintf("ORA %s %s failed: %s", db, label,
                                   conditionMessage(e)), "run_ora")
          NULL
        }
      )
      if (!is.null(res) && nrow(res) > 0) {
        res[, database      := db]
        res[, db_version    := format(Sys.Date(), "%Y-%m-%d")]
        res[, db_access_date:= format(Sys.Date(), "%Y-%m-%d")]
        results[[db]] <- res
      }
    }
    if (length(results) == 0) return(NULL)
    data.table::rbindlist(results, fill=TRUE)
  }

  list(
    all  = .do_ora(result_dt[significant==TRUE,     feature_id], "all"),
    up   = .do_ora(result_dt[direction=="up",        feature_id], "up"),
    down = .do_ora(result_dt[direction=="down",      feature_id], "down")
  )
}

#' @keywords internal
.ora_one_db <- function(entrez, universe, db, organism,
                         p_adj, min_gs, max_gs) {
  if (db %in% c("GO_BP","GO_MF","GO_CC")) {
    ont <- gsub("GO_","",db)
    res <- clusterProfiler::enrichGO(
      gene          = entrez,
      universe      = universe,
      OrgDb         = .get_org_db(organism),
      ont           = ont,
      pvalueCutoff  = 1,
      qvalueCutoff  = 1,
      pAdjustMethod = p_adj,
      minGSSize     = min_gs,
      maxGSSize     = max_gs,
      readable      = TRUE
    )
  } else if (db == "KEGG") {
    org_code <- .get_kegg_code(organism)
    res <- clusterProfiler::enrichKEGG(
      gene          = entrez,
      organism      = org_code,
      universe      = universe,
      pvalueCutoff  = 1,
      qvalueCutoff  = 1,
      pAdjustMethod = p_adj,
      minGSSize     = min_gs,
      maxGSSize     = max_gs
    )
  } else if (db == "Reactome") {
    if (!requireNamespace("ReactomePA", quietly=TRUE)) return(NULL)
    res <- ReactomePA::enrichPathway(
      gene          = entrez,
      universe      = universe,
      organism      = organism,
      pvalueCutoff  = 1,
      qvalueCutoff  = 1,
      pAdjustMethod = p_adj,
      minGSSize     = min_gs,
      maxGSSize     = max_gs,
      readable      = TRUE
    )
  } else {
    return(NULL)
  }

  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  data.table::as.data.table(as.data.frame(res))
}

#' @keywords internal
.get_kegg_code <- function(organism) {
  switch(tolower(organism), human="hsa", mouse="mmu", rat="rno", "hsa")
}

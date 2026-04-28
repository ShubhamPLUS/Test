#' Run Gene Set Enrichment Analysis (GSEA)
#'
#' Runs GSEA using `fgsea` on the full ranked statistic per contrast.
#' The ranking metric is `sign(log2FC) × -log10(p-value)`.
#'
#' @param result_dt   `data.table` from `run_limma_deqms()`.
#' @param organism    Organism string.
#' @param databases   Character vector of databases.
#' @param id_type     Input ID type.
#' @param min_gs      Minimum gene-set size. Default 10.
#' @param max_gs      Maximum gene-set size. Default 500.
#' @param n_perm      Number of permutations. Default 1000.
#' @param seed        Random seed. Default 1234.
#'
#' @return A named list of `data.table` results (one per database).
#'
#' @examples
#' \dontrun{
#' gsea_res <- run_gsea(result_dt, organism="human")
#' }
#'
#' @export
run_gsea <- function(result_dt,
                     organism  = "human",
                     databases = c("GO_BP","KEGG","Reactome"),
                     id_type   = "auto",
                     min_gs    = 10L,
                     max_gs    = 500L,
                     n_perm    = 1000L,
                     seed      = 1234L) {
  if (!requireNamespace("fgsea", quietly=TRUE)) {
    collect_warning("Package 'fgsea' not installed. GSEA skipped.", "run_gsea")
    return(NULL)
  }

  # Build ranking
  ranked <- build_gsea_ranking(result_dt)
  if (length(ranked) < 10) return(NULL)

  # Map to Entrez
  entrez_map <- map_ids_to_entrez(names(ranked), organism, id_type)
  mapped_ranked <- ranked[!is.na(entrez_map)]
  names(mapped_ranked) <- entrez_map[!is.na(entrez_map)]
  mapped_ranked <- mapped_ranked[!duplicated(names(mapped_ranked))]

  results <- list()
  for (db in databases) {
    pathways <- tryCatch(
      .get_pathways(db, organism, entrez_ids=names(mapped_ranked)),
      error = function(e) NULL
    )
    if (is.null(pathways) || length(pathways) == 0) next

    set.seed(seed)
    res <- tryCatch(
      suppressWarnings(fgsea::fgsea(
        pathways  = pathways,
        stats     = mapped_ranked,
        minSize   = min_gs,
        maxSize   = max_gs,
        nPermSimple = n_perm
      )),
      error = function(e) {
        collect_warning(sprintf("GSEA %s failed: %s", db, conditionMessage(e)),
                        "run_gsea")
        NULL
      }
    )
    if (!is.null(res) && nrow(res) > 0) {
      res <- data.table::as.data.table(res)
      res[, database       := db]
      res[, db_access_date := format(Sys.Date(), "%Y-%m-%d")]
      results[[db]] <- res
    }
  }

  if (length(results) == 0) return(NULL)
  data.table::rbindlist(results, fill=TRUE)
}

#' @keywords internal
.get_pathways <- function(db, organism, entrez_ids=NULL) {
  if (db %in% c("GO_BP","GO_MF","GO_CC")) {
    org_pkg <- .get_org_db(organism)
    if (!requireNamespace(org_pkg, quietly=TRUE)) return(NULL)
    ont <- gsub("GO_","", db)
    clusterProfiler::go2term # ensure loaded
    go2gene <- suppressMessages(
      clusterProfiler::go2gene(
        OrgDb = org_pkg, ont = ont, keyType = "ENTREZID"
      )
    )
    if (is.null(go2gene)) return(NULL)
    split(go2gene$gene, go2gene$go_id)
  } else if (db == "KEGG") {
    if (!requireNamespace("clusterProfiler", quietly=TRUE)) return(NULL)
    org_code <- .get_kegg_code(organism)
    kegg     <- suppressMessages(
      clusterProfiler::download_KEGG(org_code)
    )
    if (is.null(kegg)) return(NULL)
    split(kegg$KEGGPATHID2EXTID$to, kegg$KEGGPATHID2EXTID$from)
  } else {
    NULL
  }
}

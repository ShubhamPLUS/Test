#' Query STRINGdb for a protein interaction network
#'
#' Retrieves interactions from STRINGdb v12 for the set of significant
#' proteins, applies a score threshold, and returns node and edge tables
#' plus a `ggraph` plot.
#'
#' @param result_dt      `data.table` from `run_limma_deqms()`.
#' @param organism       Organism string. Default `"human"`.
#' @param score_threshold Combined score threshold (0–1000). Default 700.
#' @param max_nodes      Maximum nodes to include in the plot. Default 200.
#' @param output_dir     Directory to write network files. If `NULL`, files
#'   are not written.
#'
#' @return A named list: `nodes`, `edges`, `graph` (igraph object),
#'   `plot` (ggraph plot), `cytoscape_nodes`, `cytoscape_edges`.
#'
#' @examples
#' \dontrun{
#' net <- run_string_network(result_dt, organism="human")
#' }
#'
#' @export
run_string_network <- function(result_dt,
                                organism        = "human",
                                score_threshold = 700L,
                                max_nodes       = 200L,
                                output_dir      = NULL) {
  if (!requireNamespace("STRINGdb", quietly=TRUE)) {
    collect_warning("Package 'STRINGdb' not installed. Network analysis skipped.",
                    "run_string_network")
    return(NULL)
  }
  if (!requireNamespace("igraph", quietly=TRUE)) {
    collect_warning("Package 'igraph' not installed. Network analysis skipped.",
                    "run_string_network")
    return(NULL)
  }

  sig_dt <- result_dt[significant==TRUE]
  if (nrow(sig_dt) == 0) {
    logger::log_info("run_string_network: no significant proteins")
    return(NULL)
  }

  # Get gene symbols for STRING query
  id_col <- if ("gene_symbol" %in% names(sig_dt) &&
                any(!is.na(sig_dt$gene_symbol))) "gene_symbol" else "feature_id"
  query_ids <- unique(na.omit(sig_dt[[id_col]]))
  if (length(query_ids) > max_nodes) {
    # Take top max_nodes by significance
    query_ids <- sig_dt[order(adj.P.Val)][seq_len(max_nodes)][[id_col]]
  }

  # STRING species code
  species_code <- switch(tolower(organism),
    human=9606L, mouse=10090L, rat=10116L, yeast=4932L, 9606L)

  db <- tryCatch({
    STRINGdb::STRINGdb$new(
      version          = "12",
      species          = species_code,
      score_threshold  = score_threshold,
      network_type     = "full",
      input_directory  = tempdir()
    )
  }, error = function(e) {
    collect_warning(sprintf("STRINGdb init failed: %s. Check network access.",
                            conditionMessage(e)), "run_string_network")
    return(NULL)
  })
  if (is.null(db)) return(NULL)

  # Map and retrieve network
  ids_df <- data.frame(gene=query_ids)
  mapped <- tryCatch(db$map(ids_df, "gene", removeUnmappedRows=TRUE),
                     error=function(e) NULL)
  if (is.null(mapped) || nrow(mapped)==0) return(NULL)

  string_ids <- mapped$STRING_id
  net <- tryCatch(
    db$get_subnetwork(string_ids),
    error = function(e) NULL
  )
  if (is.null(net)) return(NULL)

  # ── Nodes ──────────────────────────────────────────────────────���───────────
  nodes_dt <- data.table::as.data.table(igraph::as_data_frame(net, what="vertices"))
  nodes_dt[, degree      := igraph::degree(net)]
  nodes_dt[, betweenness := igraph::betweenness(net, normalized=TRUE)]
  nodes_dt[, eigenvector := igraph::eigen_centrality(net)$vector]

  # Merge DE stats
  merge_cols <- intersect(c("feature_id","gene_symbol","log2FC","adj.P.Val","direction"),
                           names(sig_dt))
  nodes_dt <- merge(nodes_dt, sig_dt[, ..merge_cols],
                     by.x="name", by.y=id_col, all.x=TRUE)

  # ── Edges ───────────────────────────────────────────────────────────��──────
  edges_dt <- data.table::as.data.table(igraph::as_data_frame(net, what="edges"))
  edges_dt[, combined_score := score_threshold]

  # ── Cytoscape export ───────────────────────────────────────────────────────
  cyto_nodes <- nodes_dt[, .(id=name, label=name, log2FC, adj.P.Val,
                               degree, betweenness)]
  cyto_edges <- edges_dt[, .(source=from, target=to)]

  # ── ggraph plot ────────────────────────────────────────────────────────────
  plot_obj <- NULL
  if (requireNamespace("ggraph", quietly=TRUE) && igraph::vcount(net) >= 3) {
    set.seed(1234)
    tryCatch({
      plot_obj <- ggraph::ggraph(net, layout="fr") +
        ggraph::geom_edge_link(colour="grey70", alpha=0.5, linewidth=0.3) +
        ggraph::geom_node_point(
          ggplot2::aes(size=degree,
                        colour=ifelse(name %in% sig_dt[direction=="up", feature_id],
                                      "up", ifelse(
                                        name %in% sig_dt[direction=="down", feature_id],
                                        "down","ns"))),
          alpha=0.9
        ) +
        ggplot2::scale_colour_manual(
          values=c(up="#D55E00",down="#0072B2",ns="grey60"),
          name="Direction"
        ) +
        ggraph::geom_node_text(ggplot2::aes(label=name), size=2.5,
                                repel=TRUE, max.overlaps=15) +
        theme_proteoforge() +
        ggplot2::labs(title=sprintf("STRING network (score≥%d, n=%d)",
                                     score_threshold, igraph::vcount(net)))
    }, error=function(e) NULL)
  }

  # Write files
  if (!is.null(output_dir)) {
    fs::dir_create(output_dir)
    data.table::fwrite(nodes_dt, fs::path(output_dir,"nodes.tsv"),  sep="\t")
    data.table::fwrite(edges_dt, fs::path(output_dir,"edges.tsv"),  sep="\t")
    data.table::fwrite(cyto_nodes, fs::path(output_dir,"cytoscape_nodes.csv"))
    data.table::fwrite(cyto_edges, fs::path(output_dir,"cytoscape_edges.csv"))
    tryCatch(
      igraph::write_graph(net, fs::path(output_dir,"network.graphml"),
                          format="graphml"),
      error=function(e) NULL
    )
    logger::log_info("Network files written to {output_dir}")
  }

  list(
    nodes           = nodes_dt,
    edges           = edges_dt,
    graph           = net,
    plot            = plot_obj,
    cytoscape_nodes = cyto_nodes,
    cytoscape_edges = cyto_edges
  )
}

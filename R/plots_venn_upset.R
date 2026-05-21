#' Venn / UpSet diagram for contrast overlaps
#'
#' Supports 2–7 sets. Uses `ggVennDiagram` for 2–4 sets and falls back to
#' `ComplexUpset` for 5+ sets (or when explicitly requested).
#'
#' @param results_list Named list of contrast result `data.table`s from
#'   `run_limma_deqms()`.
#' @param mode   One of `"all"`, `"up"`, `"down"`.
#' @param type   One of `"venn"` (auto), `"upset"`.
#' @param title  Plot title.
#'
#' @return A ggplot2-compatible plot object.
#'
#' @examples
#' \dontrun{
#' p <- plot_overlap(list(A = dt1, B = dt2), mode = "all")
#' }
#'
#' @export
plot_overlap <- function(results_list, mode = "all", type = "auto",
                         title = NULL) {
  n <- length(results_list)
  if (n < 2) stop("Need at least 2 contrasts for overlap analysis.", call.=FALSE)

  # Extract significant IDs
  sets <- lapply(results_list, function(dt) {
    if (mode == "up")   return(dt[direction == "up",   feature_id])
    if (mode == "down") return(dt[direction == "down", feature_id])
    dt[significant == TRUE, feature_id]
  })

  use_upset <- (type == "upset") || (n > 4)

  if (use_upset) {
    .plot_upset(sets, title %||% sprintf("UpSet — %s", mode))
  } else {
    .plot_venn(sets, title %||% sprintf("Venn — %s", mode))
  }
}

#' @keywords internal
.plot_venn <- function(sets, title) {
  if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
    stop("Package 'ggVennDiagram' required for Venn diagrams.", call.=FALSE)
  }
  p <- ggVennDiagram::ggVennDiagram(sets, label_alpha = 0) +
    ggplot2::scale_fill_gradient(low = "white", high = "#0072B2") +
    ggplot2::ggtitle(title) +
    theme_proteoforge()
  p
}

#' @keywords internal
.plot_upset <- function(sets, title) {
  if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
    stop("Package 'ComplexUpset' required for UpSet plots.", call.=FALSE)
  }
  all_ids <- unique(unlist(sets))
  dt <- data.table::data.table(feature_id = all_ids)
  for (nm in names(sets)) {
    dt[[nm]] <- dt$feature_id %in% sets[[nm]]
  }
  set_cols <- setdiff(names(dt), "feature_id")
  ComplexUpset::upset(as.data.frame(dt), set_cols,
                       name = title, width_ratio = 0.15) +
    ggplot2::ggtitle(title)
}

#' Extract overlap intersection table
#'
#' @param results_list Named list of contrast result data.tables.
#' @param mode         `"all"`, `"up"`, or `"down"`.
#'
#' @return A `data.table` with columns `intersection_id`, `set_membership`,
#'   `feature_ids`, `gene_symbols`, `count`.
#'
#' @examples
#' dt1 <- data.table::data.table(feature_id=c("P1","P2","P3"),
#'   gene_symbol=c("G1","G2","G3"), significant=TRUE, direction="up")
#' dt2 <- data.table::data.table(feature_id=c("P2","P3","P4"),
#'   gene_symbol=c("G2","G3","G4"), significant=TRUE, direction="up")
#' get_intersection_table(list(A=dt1, B=dt2))
#'
#' @export
get_intersection_table <- function(results_list, mode = "all") {
  sets <- lapply(results_list, function(dt) {
    if (mode == "up")   return(dt[direction == "up",   feature_id])
    if (mode == "down") return(dt[direction == "down", feature_id])
    dt[significant == TRUE, feature_id]
  })

  all_ids  <- unique(unlist(sets))
  set_nms  <- names(sets)

  # Build membership matrix
  mem <- vapply(set_nms, function(s) all_ids %in% sets[[s]], logical(length(all_ids)))
  rownames(mem) <- all_ids

  # Group by membership pattern
  patterns <- apply(mem, 1, function(r) paste(set_nms[r], collapse="&"))

  dt_out <- data.table::data.table(
    feature_id    = all_ids,
    set_membership= patterns
  )

  # Attach gene symbols from first result list element
  first_dt  <- results_list[[1]]
  gs_col    <- if ("gene_symbol" %in% names(first_dt)) "gene_symbol" else "feature_id"
  all_genes <- rbind(rbindlist(lapply(results_list, function(d)
    d[, .(feature_id, gene_symbol = if (gs_col %in% names(d)) d[[gs_col]] else feature_id)]),
    fill = TRUE))
  all_genes <- unique(all_genes)
  dt_out    <- merge(dt_out, all_genes, by = "feature_id", all.x = TRUE)

  # Summarise by intersection
  summary_dt <- dt_out[, .(
    feature_ids  = paste(feature_id,  collapse=";"),
    gene_symbols = paste(gene_symbol, collapse=";"),
    count        = .N
  ), by = set_membership]
  summary_dt[, intersection_id := paste0("INT_", .I)]
  data.table::setcolorder(summary_dt, c("intersection_id","set_membership",
                                         "feature_ids","gene_symbols","count"))
  data.table::setorder(summary_dt, -count)
  summary_dt
}

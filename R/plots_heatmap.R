#' Heatmap of top N differentially abundant features
#'
#' Uses `pheatmap` for a clustered heatmap of Z-scored log2 intensities.
#'
#' @param pfd        A `ProteoForgeData` object.
#' @param result_dt  `data.table` from `run_limma_deqms()`.
#' @param top_n      Number of top features to include. Default 50.
#' @param label_col  Feature metadata column for row labels. Default `"gene_symbol"`.
#' @param title      Plot title.
#'
#' @return A `pheatmap` object.
#'
#' @examples
#' \dontrun{
#' p <- plot_heatmap(pfd, result_dt, top_n = 30)
#' }
#'
#' @export
plot_heatmap <- function(pfd, result_dt,
                         top_n     = 50L,
                         label_col = "gene_symbol",
                         title     = "Top DE features") {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' required.", call. = FALSE)
  }

  mat  <- pfd@raw_matrix
  smd  <- pfd@sample_metadata
  fmd  <- pfd@feature_metadata

  # Select top features
  data.table::setorder(result_dt, adj.P.Val, -abs(log2FC))
  top_ids <- result_dt[seq_len(min(top_n, nrow(result_dt))), feature_id]
  top_ids <- intersect(top_ids, rownames(mat))
  if (length(top_ids) == 0) {
    warning("No features to plot in heatmap.", call. = FALSE)
    return(invisible(NULL))
  }

  sub_mat <- mat[top_ids, , drop = FALSE]

  # Z-score rows
  sub_mat_z <- t(scale(t(sub_mat)))

  # Row labels
  if (!is.null(fmd) && label_col %in% names(fmd)) {
    lab_idx  <- match(top_ids, fmd$feature_id)
    row_labs <- ifelse(is.na(lab_idx), top_ids, fmd[[label_col]][lab_idx])
    row_labs[is.na(row_labs)] <- top_ids[is.na(row_labs)]
  } else {
    row_labs <- top_ids
  }
  rownames(sub_mat_z) <- row_labs

  # Column annotation
  ann_col <- as.data.frame(smd[, .(condition)], row.names = smd$sample_id)
  n_conds <- length(unique(smd$condition))
  ann_colours <- list(condition = stats::setNames(pf_palette(n_conds),
                                                   unique(smd$condition)))

  pheatmap::pheatmap(
    sub_mat_z,
    color            = grDevices::colorRampPalette(
                         c("#0072B2","white","#D55E00"))(100),
    annotation_col   = ann_col,
    annotation_colors= ann_colours,
    cluster_rows     = TRUE,
    cluster_cols     = TRUE,
    show_rownames    = TRUE,
    show_colnames    = TRUE,
    fontsize_row     = max(5, 9 - nrow(sub_mat_z) %/% 10),
    main             = title,
    silent           = TRUE
  )
}

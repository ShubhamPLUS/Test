#' Volcano plot for differential abundance results
#'
#' @param result_dt   `data.table` from `run_limma_deqms()`.
#' @param title       Plot title.
#' @param alpha       Adjusted p-value threshold (horizontal dashed line).
#' @param lfc_cutoff  Absolute log2FC threshold (vertical dashed lines).
#' @param label_top_n Number of top significant features to label.
#' @param label_col   Column to use for labels. Default `"gene_symbol"`.
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' dt <- data.table::data.table(
#'   feature_id  = paste0("P", 1:100),
#'   gene_symbol = paste0("G", 1:100),
#'   log2FC      = rnorm(100, 0, 2),
#'   adj.P.Val   = runif(100, 0, 1),
#'   P.Value     = runif(100, 0, 1),
#'   significant = FALSE,
#'   direction   = "ns"
#' )
#' plot_volcano(dt)
#'
#' @export
plot_volcano <- function(result_dt,
                         title      = "Volcano Plot",
                         alpha      = 0.05,
                         lfc_cutoff = 1,
                         label_top_n = 15,
                         label_col   = "gene_symbol") {

  dt <- data.table::copy(result_dt)
  dt[, neg_log10_p := -log10(pmax(adj.P.Val, 1e-300))]

  if (!label_col %in% names(dt)) label_col <- "feature_id"

  # Assign colour category
  dt[, colour_group := "ns"]
  dt[direction == "up",   colour_group := "up"]
  dt[direction == "down", colour_group := "down"]

  cols <- c(up = "#D55E00", down = "#0072B2", ns = "grey70")

  # Select top features to label
  dt_sig  <- dt[significant == TRUE]
  data.table::setorder(dt_sig, adj.P.Val, -abs(log2FC))
  top_ids <- dt_sig[seq_len(min(label_top_n, nrow(dt_sig))), feature_id]
  dt[, label := ifelse(feature_id %in% top_ids, get(label_col), "")]

  ggplot2::ggplot(dt, ggplot2::aes(x = log2FC, y = neg_log10_p,
                                    colour = colour_group)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.2) +
    ggplot2::scale_colour_manual(
      values = cols,
      labels = c(up   = sprintf("Up (n=%d)",   sum(dt$direction=="up")),
                 down = sprintf("Down (n=%d)", sum(dt$direction=="down")),
                 ns   = sprintf("NS (n=%d)",   sum(dt$direction=="ns"))),
      name = NULL
    ) +
    ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed",
                         colour = "grey40", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
                         linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    ggrepel::geom_text_repel(
      data = dt[label != ""],
      ggplot2::aes(label = label),
      size = 2.8, max.overlaps = 30, seed = 1234
    ) +
    ggplot2::labs(
      title = title,
      x     = expression(log[2]~"Fold Change"),
      y     = expression(-log[10]~"(adj. p-value)")
    ) +
    theme_proteoforge()
}

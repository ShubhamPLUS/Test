#' Dotplot for ORA results
#'
#' @param ora_dt   `data.table` from `run_ora()`.
#' @param top_n    Top N terms to display. Default 20.
#' @param title    Plot title.
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' \dontrun{
#' p <- plot_ora_dotplot(ora_dt, top_n=20)
#' }
#'
#' @export
plot_ora_dotplot <- function(ora_dt, top_n=20L, title="ORA Results") {
  if (is.null(ora_dt) || nrow(ora_dt)==0) return(NULL)

  dt <- data.table::copy(ora_dt)
  p_col <- intersect(c("p.adjust","qvalue","pvalue"), names(dt))[1]
  if (is.na(p_col)) return(NULL)

  dt[, neg_log_p := -log10(pmax(as.numeric(get(p_col)), 1e-300))]
  dt <- dt[order(-neg_log_p)][seq_len(min(top_n, nrow(dt)))]

  desc_col <- intersect(c("Description","pathway","ID"), names(dt))[1]
  count_col <- intersect(c("Count","overlap","size"), names(dt))[1]

  if (is.na(desc_col)) return(NULL)

  ggplot2::ggplot(dt,
    ggplot2::aes(x=neg_log_p,
                  y=stats::reorder(get(desc_col), neg_log_p),
                  size=if (!is.na(count_col)) as.numeric(get(count_col)) else 3,
                  colour=neg_log_p)) +
    ggplot2::geom_point(alpha=0.8) +
    ggplot2::scale_colour_gradient(low="#56B4E9", high="#D55E00") +
    ggplot2::labs(x=expression(-log[10]~"(adj. p-value)"),
                   y=NULL, title=title,
                   size="Gene count", colour=expression(-log[10]~p)) +
    theme_proteoforge()
}

#' Ridge plot for GSEA results
#'
#' @param gsea_dt   `data.table` from `run_gsea()`.
#' @param top_n     Top N pathways. Default 10.
#' @param title     Plot title.
#'
#' @return A `ggplot2` object or `NULL`.
#'
#' @examples
#' \dontrun{
#' p <- plot_gsea_ridge(gsea_dt)
#' }
#'
#' @export
plot_gsea_ridge <- function(gsea_dt, top_n=10L, title="GSEA Top Pathways") {
  if (is.null(gsea_dt) || nrow(gsea_dt)==0) return(NULL)
  if (!requireNamespace("ggridges", quietly=TRUE)) {
    collect_warning("Package 'ggridges' not installed. GSEA ridge plot skipped.",
                    "plot_gsea_ridge")
    return(NULL)
  }

  dt <- data.table::copy(gsea_dt)
  p_col <- intersect(c("padj","pval","p.adjust"), names(dt))[1]
  if (is.na(p_col)) return(NULL)

  dt[, signed_score := ifelse(NES > 0, -log10(pmax(as.numeric(get(p_col)),1e-300)),
                               log10(pmax(as.numeric(get(p_col)),1e-300)))]
  dt <- dt[order(-abs(signed_score))][seq_len(min(top_n, nrow(dt)))]

  desc_col <- intersect(c("pathway","Description","ID"), names(dt))[1]
  if (is.na(desc_col)) return(NULL)

  ggplot2::ggplot(dt,
    ggplot2::aes(x=signed_score,
                  y=stats::reorder(get(desc_col), signed_score),
                  fill=NES > 0)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values=c("TRUE"="#D55E00","FALSE"="#0072B2"),
                                labels=c("TRUE"="Enriched","FALSE"="Depleted"),
                                name=NULL) +
    ggplot2::labs(x="Signed –log10(adj.p)", y=NULL, title=title) +
    theme_proteoforge()
}

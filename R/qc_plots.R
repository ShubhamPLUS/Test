#' Generate QC plots from a ProteoForgeData object
#'
#' Creates the full suite of QC visualisations and saves them using
#' `save_plot()`. All plots are returned in a named list.
#'
#' @param pfd         A `ProteoForgeData` object.
#' @param qc_metrics  Named list from `compute_qc_metrics()`.
#' @param plot_dir    Directory to save plots.
#' @param config      Named list following the ProteoForge config schema.
#'
#' @return Named list of ggplot2 objects.
#'
#' @examples
#' \dontrun{
#' qf  <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
#'                    package="proteoforge")
#' mf  <- system.file("extdata","generic_lfq_small","sample_metadata.tsv",
#'                    package="proteoforge")
#' pfd <- import_generic(qf, mf)
#' pfd <- log2_transform(pfd)
#' pfd <- normalise_matrix(pfd)
#' qc  <- compute_qc_metrics(pfd)
#' plots <- generate_qc_plots(pfd, qc, tempdir())
#' }
#'
#' @export
generate_qc_plots <- function(pfd, qc_metrics, plot_dir,
                               config = list()) {

  stopifnot(methods::is(pfd, "ProteoForgeData"))
  fs::dir_create(plot_dir)

  project  <- config$project$name       %||% "proteoforge"
  formats  <- config$plots$formats      %||% c("pdf","png")
  dpi      <- config$plots$dpi          %||% 300
  phash    <- pfd@parameter_hash

  plots <- list()

  # ── 1. ID counts per sample ───────────────────────────────────────────────
  plots$id_counts <- .plot_id_counts(qc_metrics$id_counts) +
    theme_proteoforge() +
    ggplot2::labs(title = "Detected features per sample")

  save_plot(plots$id_counts, plot_dir,
            safe_filename(project,"QC","protein","all","id_counts","",
                          Sys.time()) |> fs::path_ext_remove(),
            formats = formats, dpi = dpi, module = "QC", plot_type = "id_counts",
            params_hash = phash)

  # ── 2. Sample intensity distribution (boxplot) ────────────────────────────
  plots$intensity_box <- .plot_intensity_boxplot(pfd) +
    theme_proteoforge() +
    ggplot2::labs(title = "Intensity distribution per sample (log2)")

  save_plot(plots$intensity_box, plot_dir,
            safe_filename(project,"QC","protein","all","intensity_boxplot","",
                          Sys.time()) |> fs::path_ext_remove(),
            formats = formats, dpi = dpi, module = "QC",
            plot_type = "intensity_boxplot", params_hash = phash)

  # ── 3. PCA ─────────────────────────────────────────────────────────────────
  if (!is.null(qc_metrics$pca)) {
    plots$pca <- .plot_pca(qc_metrics$pca, pfd@sample_metadata) +
      theme_proteoforge() +
      ggplot2::labs(title = "PCA — samples coloured by condition")

    save_plot(plots$pca, plot_dir,
              safe_filename(project,"QC","protein","all","pca","",
                            Sys.time()) |> fs::path_ext_remove(),
              formats = formats, dpi = dpi, module = "QC", plot_type = "pca",
              params_hash = phash)
  }

  # ── 4. Correlation heatmap ────────────────────────────────────────────────
  if (!is.null(qc_metrics$correlation)) {
    plots$correlation <- .plot_correlation_heatmap(qc_metrics$correlation,
                                                    pfd@sample_metadata) +
      ggplot2::labs(title = "Sample Pearson correlation (log2 intensity)")

    save_plot(plots$correlation, plot_dir,
              safe_filename(project,"QC","protein","all","correlation","",
                            Sys.time()) |> fs::path_ext_remove(),
              formats = formats, dpi = dpi, module = "QC",
              plot_type = "correlation_heatmap", params_hash = phash)
  }

  # ── 5. Missingness heatmap ────────────────────────────────────────────────
  plots$missingness <- .plot_missingness(pfd) +
    theme_proteoforge() +
    ggplot2::labs(title = "Missingness by sample and feature")

  save_plot(plots$missingness, plot_dir,
            safe_filename(project,"QC","protein","all","missingness","",
                          Sys.time()) |> fs::path_ext_remove(),
            formats = formats, dpi = dpi, module = "QC",
            plot_type = "missingness", params_hash = phash)

  invisible(plots)
}

# ── Private plot helpers ──────────────────────────────────────────────────────

#' @keywords internal
.plot_id_counts <- function(id_counts) {
  ggplot2::ggplot(id_counts,
                  ggplot2::aes(x = reorder(sample_id, n_detected),
                               y = n_detected,
                               fill = condition)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    scale_colour_pf() +
    ggplot2::scale_fill_manual(values = pf_palette(
      length(unique(id_counts$condition)))) +
    ggplot2::labs(x = NULL, y = "# Detected features", fill = "Condition")
}

#' @keywords internal
.plot_intensity_boxplot <- function(pfd) {
  mat  <- pfd@raw_matrix
  smd  <- pfd@sample_metadata

  long <- data.table::as.data.table(
    reshape2::melt(mat, varnames = c("feature_id","sample_id"),
                   value.name = "intensity")
  )
  long <- merge(long, smd[, .(sample_id, condition)], by="sample_id",
                all.x=TRUE)
  long <- long[!is.na(intensity)]

  ggplot2::ggplot(long, ggplot2::aes(x = sample_id, y = intensity,
                                      fill = condition)) +
    ggplot2::geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
    ggplot2::scale_fill_manual(values = pf_palette(
      length(unique(long$condition)))) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "log2 Intensity", fill = "Condition")
}

#' @keywords internal
.plot_pca <- function(pca_result, smd) {
  pca_dt  <- pca_result$pca_dt
  var_exp <- pca_result$var_exp

  ggplot2::ggplot(pca_dt,
                  ggplot2::aes(x = PC1, y = PC2, colour = condition,
                               label = sample_id)) +
    ggplot2::geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
    ggplot2::scale_colour_manual(values = pf_palette(
      length(unique(pca_dt$condition)))) +
    ggplot2::labs(
      x      = sprintf("PC1 (%.1f%%)", var_exp[1]),
      y      = sprintf("PC2 (%.1f%%)", var_exp[2]),
      colour = "Condition"
    )
}

#' @keywords internal
.plot_correlation_heatmap <- function(cor_mat, smd) {
  # Convert to long format for ggplot2
  cor_long <- data.table::as.data.table(
    reshape2::melt(cor_mat, varnames = c("Sample1","Sample2"),
                   value.name = "Pearson_r")
  )

  ggplot2::ggplot(cor_long,
                  ggplot2::aes(x = Sample1, y = Sample2, fill = Pearson_r)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#0072B2", mid = "white",
                                   high = "#D55E00", midpoint = 0.95,
                                   limits = c(NA, 1)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1,
                                                        size=7),
                   axis.text.y = ggplot2::element_text(size=7)) +
    ggplot2::labs(x=NULL, y=NULL, fill="Pearson r") +
    ggplot2::coord_fixed()
}

#' @keywords internal
.plot_missingness <- function(pfd) {
  mat  <- pfd@raw_matrix
  smd  <- pfd@sample_metadata

  miss_long <- data.table::as.data.table(
    reshape2::melt(is.na(mat), varnames = c("feature_id","sample_id"),
                   value.name = "is_missing")
  )
  miss_long <- merge(miss_long, smd[, .(sample_id, condition)],
                     by = "sample_id", all.x = TRUE)

  # Summarise per sample
  sample_miss <- miss_long[, .(pct_missing = 100*mean(is_missing)), by=sample_id]
  sample_miss <- merge(sample_miss, smd[, .(sample_id, condition)],
                       by="sample_id", all.x=TRUE)

  ggplot2::ggplot(sample_miss,
                  ggplot2::aes(x = reorder(sample_id, pct_missing),
                               y = pct_missing, fill = condition)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = pf_palette(
      length(unique(sample_miss$condition)))) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "% Missing values", fill = "Condition") +
    ggplot2::geom_hline(yintercept = 30, linetype = "dashed",
                         colour = "red", linewidth = 0.5)
}

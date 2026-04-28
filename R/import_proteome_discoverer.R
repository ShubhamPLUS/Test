#' Import Proteome Discoverer protein abundance export
#'
#' Reads PD proteins TSV/Excel export, applies master protein and confidence
#' filters, and extracts `Abundance: F<n>: Sample, ...` columns via regex.
#' Provides a manual column mapping fallback for non-standard templates.
#'
#' @param quant_file       Path to PD proteins export (TSV or Excel).
#' @param metadata_file    Path to sample metadata.
#' @param config           ProteoForge config list.
#' @param abundance_pattern Regex to identify abundance columns. Default
#'   matches `Abundance: F\\d+:`.
#' @param sample_name_fn   Optional function `f(col_name) -> sample_id` for
#'   non-standard column name mapping.
#'
#' @return A `ProteoForgeData` object.
#'
#' @examples
#' qf <- system.file("extdata","pd_example","proteins.tsv",
#'                   package="proteoforge")
#' mf <- system.file("extdata","pd_example","sample_metadata.tsv",
#'                   package="proteoforge")
#' pfd <- import_pd(qf, mf)
#'
#' @export
import_pd <- function(quant_file, metadata_file, config = list(),
                       abundance_pattern = "^Abundance: F\\d+:",
                       sample_name_fn    = NULL) {
  stopifnot(file.exists(quant_file), file.exists(metadata_file))
  logger::log_info("import_pd: {quant_file}")

  smd <- .read_metadata(metadata_file)
  data.table::setDT(smd)
  cfg <- config$filtering %||% list()

  dt <- .read_delim_auto(quant_file)

  # ── Filters ────────────────────────────────────────────────────────────────
  if (isTRUE(cfg$require_master_protein %||% FALSE) &&
      "Master" %in% names(dt)) {
    dt <- dt[Master == "IsMasterProtein"]
  }

  if ("Protein FDR Confidence: Combined" %in% names(dt)) {
    dt <- dt[`Protein FDR Confidence: Combined` == "High" |
               is.na(`Protein FDR Confidence: Combined`)]
  }

  if ("Contaminant" %in% names(dt)) {
    dt <- dt[Contaminant != TRUE | is.na(Contaminant)]
  }

  # ── Abundance columns ──────────────────────────────────────────────────────
  abund_cols <- grep(abundance_pattern, names(dt), value=TRUE)
  if (length(abund_cols) == 0) {
    stop("No abundance columns found matching '", abundance_pattern,
         "'. Provide a custom abundance_pattern or sample_name_fn.", call.=FALSE)
  }

  if (!is.null(sample_name_fn)) {
    col_sample_map <- setNames(sapply(abund_cols, sample_name_fn), abund_cols)
  } else {
    # Default: extract file index, match to metadata by position
    file_idx <- as.integer(gsub(".*F(\\d+):.*", "\\1", abund_cols))
    col_sample_map <- setNames(smd$sample_id[file_idx], abund_cols)
  }

  shared_cols <- abund_cols[col_sample_map %in% smd$sample_id]
  if (length(shared_cols) == 0) {
    stop("No abundance columns could be matched to sample metadata.", call.=FALSE)
  }

  mat <- as.matrix(dt[, shared_cols, with=FALSE])
  mode(mat) <- "numeric"
  mat[mat == 0] <- NA
  colnames(mat) <- col_sample_map[shared_cols]

  id_col <- intersect(c("Accession","Protein Accession","Master Protein Accessions"),
                       names(dt))[1]
  if (is.na(id_col)) id_col <- names(dt)[1]
  rownames(mat) <- as.character(dt[[id_col]])

  meta_cols <- intersect(c("Accession","Description","Gene Symbol",
                             "# Peptides","# Unique Peptides",
                             "Coverage [%]","Master"),
                          names(dt))
  fmd <- dt[, meta_cols, with=FALSE]
  fmd[, feature_id := dt[[id_col]]]
  if ("Gene Symbol" %in% names(fmd)) data.table::setnames(fmd,"Gene Symbol","gene_symbol")
  if ("# Unique Peptides" %in% names(fmd)) data.table::setnames(fmd,"# Unique Peptides","n_unique_peptides")

  smd2 <- smd[sample_id %in% colnames(mat)]
  long <- data.table::melt(
    data.table::data.table(feature_id=rownames(mat), mat),
    id.vars="feature_id", variable.name="sample_id",
    value.name="intensity", variable.factor=FALSE)

  logger::log_info("import_pd: {nrow(mat)} proteins × {ncol(mat)} samples")

  ProteoForgeData(
    raw_long        = long,
    raw_matrix      = mat,
    feature_metadata= fmd,
    sample_metadata = smd2,
    analysis_level  = config$input$analysis_level %||% "protein",
    quant_type      = "intensity",
    source_software = "pd",
    source_files    = c(quant_file, metadata_file),
    parameters      = list(abundance_pattern=abundance_pattern)
  )
}

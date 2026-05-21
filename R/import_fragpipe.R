#' Import FragPipe / MSFragger / IonQuant combined_protein.tsv
#'
#' Reads `combined_protein.tsv`, applies probability filters, and extracts
#' MaxLFQ Intensity columns via regex.
#'
#' @param quant_file    Path to `combined_protein.tsv`.
#' @param metadata_file Path to sample metadata.
#' @param config        ProteoForge config list.
#'
#' @return A `ProteoForgeData` object.
#'
#' @examples
#' qf <- system.file("extdata","fragpipe_example","combined_protein.tsv",
#'                   package="proteoforge")
#' mf <- system.file("extdata","fragpipe_example","sample_metadata.tsv",
#'                   package="proteoforge")
#' pfd <- import_fragpipe(qf, mf)
#'
#' @export
import_fragpipe <- function(quant_file, metadata_file, config = list()) {
  stopifnot(file.exists(quant_file), file.exists(metadata_file))
  logger::log_info("import_fragpipe: {quant_file}")

  smd <- .read_metadata(metadata_file)
  data.table::setDT(smd)
  cfg <- config$filtering %||% list()

  dt <- .read_delim_auto(quant_file)

  # ── Probability filters ────────────────────────────────────────────────────
  min_pp <- cfg$min_protein_probability %||% 0.99
  if ("Protein Probability" %in% names(dt)) {
    dt <- dt[suppressWarnings(as.numeric(`Protein Probability`)) >= min_pp |
               is.na(`Protein Probability`)]
  }

  min_pep_p <- cfg$min_top_peptide_probability %||% 0.9
  if ("Top Peptide Probability" %in% names(dt)) {
    dt <- dt[suppressWarnings(as.numeric(`Top Peptide Probability`)) >= min_pep_p |
               is.na(`Top Peptide Probability`)]
  }

  min_uniq_pep <- cfg$min_unique_peptides %||% 1L
  pep_col_name <- intersect(c("Unique Peptide Count","Unique peptides"),
                             names(dt))[1]
  if (!is.na(pep_col_name)) {
    dt <- dt[suppressWarnings(as.integer(dt[[pep_col_name]])) >= min_uniq_pep |
               is.na(dt[[pep_col_name]])]
  }

  # ── Identify MaxLFQ intensity columns ─────────────────────────────────────
  lfq_cols <- grep("MaxLFQ Intensity$", names(dt), value=TRUE)
  if (length(lfq_cols) == 0) {
    lfq_cols <- grep("[Ii]ntensity", names(dt), value=TRUE)
    lfq_cols <- setdiff(lfq_cols, c("Intensity","Total Intensity"))
  }

  # Sample name = everything before " MaxLFQ Intensity"
  col_to_sample <- gsub(" MaxLFQ Intensity$", "", lfq_cols)
  col_to_sample <- gsub(" Intensity$", "", col_to_sample)

  shared_samples <- intersect(smd$sample_id, col_to_sample)
  if (length(shared_samples) == 0) {
    stop("No sample IDs match FragPipe MaxLFQ Intensity column names.", call.=FALSE)
  }
  matched_cols <- lfq_cols[col_to_sample %in% shared_samples]

  mat <- as.matrix(dt[, matched_cols, with=FALSE])
  mode(mat) <- "numeric"
  mat[mat == 0] <- NA
  colnames(mat) <- col_to_sample[col_to_sample %in% shared_samples]

  # Protein ID column
  id_col <- intersect(c("Protein","Protein ID"), names(dt))[1]
  if (is.na(id_col)) id_col <- names(dt)[1]
  rownames(mat) <- as.character(dt[[id_col]])

  # Feature metadata
  meta_cols <- intersect(c("Protein","Protein ID","Entry Name","Gene",
                             "Protein Description","Unique Peptide Count",
                             "Razor Peptide Count"),
                          names(dt))
  fmd <- dt[, meta_cols, with=FALSE]
  data.table::setnames(fmd, id_col, "feature_id")
  if ("Gene" %in% names(fmd)) {
    data.table::setnames(fmd, "Gene", "gene_symbol")
  }
  if ("Unique Peptide Count" %in% names(fmd)) {
    data.table::setnames(fmd, "Unique Peptide Count", "n_unique_peptides")
  }

  smd2 <- smd[sample_id %in% shared_samples]
  long <- data.table::melt(
    data.table::data.table(feature_id=rownames(mat), mat),
    id.vars="feature_id", variable.name="sample_id",
    value.name="intensity", variable.factor=FALSE)

  logger::log_info("import_fragpipe: {nrow(mat)} proteins × {ncol(mat)} samples")

  ProteoForgeData(
    raw_long        = long,
    raw_matrix      = mat,
    feature_metadata= fmd,
    sample_metadata = smd2,
    analysis_level  = config$input$analysis_level %||% "protein",
    quant_type      = "maxlfq",
    source_software = "fragpipe",
    source_files    = c(quant_file, metadata_file),
    parameters      = list(quant_col="MaxLFQ Intensity",
                            min_protein_probability=min_pp)
  )
}

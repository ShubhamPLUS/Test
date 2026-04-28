#' Import DIA-NN quantification output
#'
#' Supports `report.pg_matrix.tsv`, `report.pr_matrix.tsv`,
#' `report.unique_genes_matrix.tsv`, and `report.tsv` / `report.parquet`.
#' Applies the canonical DIA-NN q-value filters.
#'
#' @param quant_file    Path to the DIA-NN output file.
#' @param metadata_file Path to sample metadata (TSV/CSV/Excel).
#' @param config        ProteoForge config list.
#'
#' @return A `ProteoForgeData` object.
#'
#' @examples
#' qf <- system.file("extdata","diann_example","report.pg_matrix.tsv",
#'                   package="proteoforge")
#' mf <- system.file("extdata","diann_example","sample_metadata.tsv",
#'                   package="proteoforge")
#' pfd <- import_diann(qf, mf)
#'
#' @export
import_diann <- function(quant_file, metadata_file, config = list()) {
  stopifnot(file.exists(quant_file), file.exists(metadata_file))
  logger::log_info("import_diann: {quant_file}")

  smd   <- .read_metadata(metadata_file)
  data.table::setDT(smd)
  assert_columns(smd, c("sample_id","condition"), "sample metadata")

  ext         <- tolower(fs::path_ext(quant_file))
  basename_lc <- tolower(fs::path_file(quant_file))

  # ── pg_matrix / pr_matrix / genes_matrix ──────────────────────────────────
  if (grepl("_matrix\\.tsv$|_matrix\\.parquet$", basename_lc)) {
    return(.import_diann_matrix(quant_file, smd, config))
  }

  # ── Full precursor-level report ────────────────────────────────────────────
  .import_diann_report(quant_file, smd, config)
}

#' @keywords internal
.import_diann_matrix <- function(path, smd, config) {
  dt <- .read_delim_auto(path)

  # Identify protein/gene annotation columns
  id_col <- intersect(c("Protein.Group","Protein.Ids","Genes"), names(dt))[1]
  if (is.na(id_col)) id_col <- names(dt)[1]

  meta_cols <- intersect(
    c("Protein.Group","Protein.Ids","Protein.Names","Genes",
      "First.Protein.Description","Protein.Description"),
    names(dt)
  )

  sample_cols <- setdiff(names(dt), meta_cols)
  # Keep only sample cols that exist in metadata
  shared <- intersect(smd$sample_id, sample_cols)
  if (length(shared) == 0) {
    stop("No sample IDs from metadata found in DIA-NN pg_matrix columns.", call.=FALSE)
  }

  mat <- as.matrix(dt[, shared, with=FALSE])
  mode(mat) <- "numeric"
  mat[mat == 0] <- NA
  rownames(mat) <- as.character(dt[[id_col]])

  fmd <- dt[, meta_cols, with=FALSE]
  data.table::setnames(fmd, id_col, "feature_id")
  if (!"gene_symbol" %in% names(fmd)) {
    g_col <- intersect(c("Genes","Protein.Names"), names(fmd))[1]
    if (!is.na(g_col)) fmd[, gene_symbol := fmd[[g_col]]]
  }

  smd2 <- smd[sample_id %in% shared]

  long <- data.table::melt(
    data.table::data.table(feature_id=rownames(mat), mat),
    id.vars="feature_id", variable.name="sample_id",
    value.name="intensity", variable.factor=FALSE)

  ProteoForgeData(
    raw_long        = long,
    raw_matrix      = mat,
    feature_metadata= fmd,
    sample_metadata = smd2,
    analysis_level  = config$input$analysis_level %||% "protein",
    quant_type      = "maxlfq",
    source_software = "diann",
    source_files    = c(path, "metadata"),
    parameters      = list(quant_file=path, zero_to_na=TRUE)
  )
}

#' @keywords internal
.import_diann_report <- function(path, smd, config) {
  qval_cutoff  <- config$filtering$qvalue_cutoff %||% 0.01
  quant_prefs  <- config$input$quant_column_preference %||%
    c("PG.MaxLFQ","Precursor.Normalised")

  # For large files, read selectively
  needed_cols <- c("Run","Protein.Group","Protein.Ids","Protein.Names","Genes",
                   "Global.Q.Value","PG.Q.Value","Lib.Q.Value","Lib.PG.Q.Value",
                   quant_prefs)

  dt <- .read_delim_auto(path)

  # Apply q-value filters
  for (qcol in intersect(c("Global.Q.Value","PG.Q.Value","Lib.Q.Value",
                             "Lib.PG.Q.Value"), names(dt))) {
    dt <- dt[is.na(get(qcol)) | get(qcol) <= qval_cutoff]
  }

  # Choose quant column
  quant_col <- Filter(function(c) c %in% names(dt), quant_prefs)[1]
  if (is.na(quant_col)) {
    stop("No recognised quant column found in DIA-NN report. ",
         "Columns present: ", paste(names(dt)[1:10], collapse=", "), call.=FALSE)
  }

  run_col <- if ("Run" %in% names(dt)) "Run" else names(dt)[1]

  # Pivot to wide
  pg_col  <- "Protein.Group"
  wide    <- data.table::dcast(dt, formula(paste(pg_col, "~", run_col)),
                                value.var = quant_col,
                                fun.aggregate = mean, na.rm=TRUE)

  # Remove placeholder rows
  wide <- wide[!is.na(get(pg_col)) & get(pg_col) != ""]

  shared    <- intersect(smd$sample_id, setdiff(names(wide), pg_col))
  mat       <- as.matrix(wide[, shared, with=FALSE])
  mode(mat) <- "numeric"
  mat[mat == 0] <- NA
  rownames(mat) <- wide[[pg_col]]

  fmd <- data.table::data.table(feature_id = wide[[pg_col]],
                                 protein_id = wide[[pg_col]])
  if ("Genes" %in% names(dt)) {
    gene_map <- unique(dt[, .(Protein.Group, gene_symbol = Genes)])
    fmd <- merge(fmd, gene_map, by.x="feature_id",
                 by.y="Protein.Group", all.x=TRUE)
  }

  smd2 <- smd[sample_id %in% shared]

  long <- data.table::melt(
    data.table::data.table(feature_id=rownames(mat), mat),
    id.vars="feature_id", variable.name="sample_id",
    value.name="intensity", variable.factor=FALSE)

  logger::log_info("import_diann: {nrow(mat)} proteins × {ncol(mat)} samples")

  ProteoForgeData(
    raw_long        = long,
    raw_matrix      = mat,
    feature_metadata= fmd,
    sample_metadata = smd2,
    analysis_level  = config$input$analysis_level %||% "protein",
    quant_type      = ifelse(quant_col == "PG.MaxLFQ", "maxlfq", "lfq"),
    source_software = "diann",
    source_files    = c(path, "metadata"),
    parameters      = list(quant_file=path, quant_col=quant_col,
                            qval_cutoff=qval_cutoff)
  )
}

#' Import Spectronaut long-format report
#'
#' Reads a Spectronaut TSV/CSV/parquet report, applies q-value filters,
#' picks the top-ranked precursor per peptide per sample, and pivots to wide.
#'
#' @param quant_file    Path to the Spectronaut report.
#' @param metadata_file Path to sample metadata.
#' @param config        ProteoForge config list.
#'
#' @return A `ProteoForgeData` object.
#'
#' @examples
#' qf <- system.file("extdata","spectronaut_example","spectronaut_report.tsv",
#'                   package="proteoforge")
#' mf <- system.file("extdata","spectronaut_example","sample_metadata.tsv",
#'                   package="proteoforge")
#' pfd <- import_spectronaut(qf, mf)
#'
#' @export
import_spectronaut <- function(quant_file, metadata_file, config = list()) {
  stopifnot(file.exists(quant_file), file.exists(metadata_file))
  logger::log_info("import_spectronaut: {quant_file}")

  smd          <- .read_metadata(metadata_file)
  data.table::setDT(smd)
  qval_cutoff  <- config$filtering$qvalue_cutoff %||% 0.01

  dt <- .read_delim_auto(quant_file)

  # ── Q-value filters ────────────────────────────────────────────────────────
  for (qcol in intersect(c("EG.Qvalue","PG.Qvalue"), names(dt))) {
    dt <- dt[is.na(get(qcol)) | get(qcol) <= qval_cutoff]
  }

  # Remove decoys
  if ("EG.IsDecoy" %in% names(dt)) {
    dt <- dt[EG.IsDecoy == FALSE | is.na(EG.IsDecoy)]
  }

  # ── Identify columns ───────────────────────────────────────────────────────
  run_col  <- intersect(c("R.FileName","R.Replicate","Sample"), names(dt))[1]
  pg_col   <- intersect(c("PG.ProteinGroups","Protein.Group"), names(dt))[1]
  quant_prefs <- config$input$quant_column_preference %||%
    c("PG.Quantity","PEP.Quantity","FG.MS2Quantity")
  quant_col <- Filter(function(c) c %in% names(dt), quant_prefs)[1]

  if (is.na(run_col))   stop("Cannot find run/filename column in Spectronaut report.", call.=FALSE)
  if (is.na(pg_col))    stop("Cannot find protein group column.", call.=FALSE)
  if (is.na(quant_col)) stop("Cannot find quant column.", call.=FALSE)

  # ── Pivot to protein × sample ──────────────────────────────────────────────
  wide <- data.table::dcast(dt,
    formula(paste(pg_col, "~", run_col)),
    value.var    = quant_col,
    fun.aggregate= function(x) mean(x, na.rm=TRUE))

  wide <- wide[!is.na(get(pg_col)) & get(pg_col) != ""]

  shared    <- intersect(smd$sample_id, setdiff(names(wide), pg_col))
  mat       <- as.matrix(wide[, shared, with=FALSE])
  mode(mat) <- "numeric"
  mat[mat == 0] <- NA
  rownames(mat) <- wide[[pg_col]]

  fmd <- data.table::data.table(feature_id = wide[[pg_col]],
                                 protein_id = wide[[pg_col]])
  if ("PG.Genes" %in% names(dt)) {
    gene_map <- unique(dt[, .(PG.ProteinGroups, gene_symbol = PG.Genes)])
    fmd <- merge(fmd, gene_map, by.x="feature_id",
                 by.y="PG.ProteinGroups", all.x=TRUE)
  }

  smd2 <- smd[sample_id %in% shared]
  long <- data.table::melt(
    data.table::data.table(feature_id=rownames(mat), mat),
    id.vars="feature_id", variable.name="sample_id",
    value.name="intensity", variable.factor=FALSE)

  logger::log_info("import_spectronaut: {nrow(mat)} proteins × {ncol(mat)} samples")

  ProteoForgeData(
    raw_long        = long,
    raw_matrix      = mat,
    feature_metadata= fmd,
    sample_metadata = smd2,
    analysis_level  = config$input$analysis_level %||% "protein",
    quant_type      = "lfq",
    source_software = "spectronaut",
    source_files    = c(quant_file, metadata_file),
    parameters      = list(quant_col=quant_col, qval_cutoff=qval_cutoff)
  )
}

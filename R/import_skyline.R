#' Import Skyline custom report export
#'
#' Skyline reports are user-defined, so column names vary widely. This
#' importer treats the file as a generic long-format report and requires
#' the user to specify column mappings (or relies on auto-detection for
#' common Skyline report templates).
#'
#' @param quant_file     Path to the Skyline custom report (TSV/CSV).
#' @param metadata_file  Path to sample metadata.
#' @param config         ProteoForge config list.
#' @param protein_col    Column name for protein/peptide ID. If `NULL`,
#'   auto-detects.
#' @param replicate_col  Column name for replicate/sample. If `NULL`,
#'   auto-detects.
#' @param quant_col      Column name for quantification value. If `NULL`,
#'   tries `"Total Area Fragment"` then `"Normalized Area"`.
#'
#' @return A `ProteoForgeData` object.
#'
#' @examples
#' \dontrun{
#' pfd <- import_skyline("skyline_report.tsv", "metadata.tsv")
#' }
#'
#' @export
import_skyline <- function(quant_file, metadata_file, config = list(),
                            protein_col   = NULL,
                            replicate_col = NULL,
                            quant_col     = NULL) {
  stopifnot(file.exists(quant_file), file.exists(metadata_file))
  logger::log_info("import_skyline: {quant_file}")

  smd <- .read_metadata(metadata_file)
  data.table::setDT(smd)

  dt  <- .read_delim_auto(quant_file)

  # ── Auto-detect columns ────────────────────────────────────────────────────
  protein_col <- protein_col %||%
    intersect(c("Protein Name","Protein","ProteinName","PeptideSequence",
                "Peptide Sequence","Modified Sequence"),
              names(dt))[1]
  if (is.na(protein_col)) protein_col <- names(dt)[1]

  replicate_col <- replicate_col %||%
    intersect(c("Replicate Name","Replicate","SampleName","File Name",
                "FileName","BioReplicate"),
              names(dt))[1]

  quant_col <- quant_col %||%
    intersect(c("Total Area Fragment","Normalized Area","Area",
                "Precursor Area","Total Area"),
              names(dt))[1]

  if (is.na(replicate_col)) {
    stop("Cannot find replicate/sample column in Skyline report. ",
         "Please specify replicate_col.", call.=FALSE)
  }
  if (is.na(quant_col)) {
    stop("Cannot find quantification column in Skyline report. ",
         "Please specify quant_col.", call.=FALSE)
  }

  collect_warning(
    sprintf("Skyline import: protein_col=%s, replicate_col=%s, quant_col=%s",
            protein_col, replicate_col, quant_col),
    "import_skyline"
  )

  # ── Pivot to wide ──────────────────────────────────────────────────────────
  wide <- data.table::dcast(dt,
    formula(paste(protein_col, "~", replicate_col)),
    value.var = quant_col,
    fun.aggregate = function(x) mean(x, na.rm=TRUE))

  shared    <- intersect(smd$sample_id, setdiff(names(wide), protein_col))
  mat       <- as.matrix(wide[, shared, with=FALSE])
  mode(mat) <- "numeric"
  mat[mat == 0] <- NA
  rownames(mat) <- as.character(wide[[protein_col]])

  fmd  <- data.table::data.table(feature_id = wide[[protein_col]],
                                   protein_id = wide[[protein_col]])
  smd2 <- smd[sample_id %in% shared]
  long <- data.table::melt(
    data.table::data.table(feature_id=rownames(mat), mat),
    id.vars="feature_id", variable.name="sample_id",
    value.name="intensity", variable.factor=FALSE)

  logger::log_info("import_skyline: {nrow(mat)} features × {ncol(mat)} samples")

  ProteoForgeData(
    raw_long        = long,
    raw_matrix      = mat,
    feature_metadata= fmd,
    sample_metadata = smd2,
    analysis_level  = config$input$analysis_level %||% "protein",
    quant_type      = "area",
    source_software = "skyline",
    source_files    = c(quant_file, metadata_file),
    parameters      = list(quant_col=quant_col, protein_col=protein_col)
  )
}

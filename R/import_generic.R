#' Import a generic feature × sample intensity matrix
#'
#' Reads a delimited file where the first column is a feature identifier and
#' remaining columns are per-sample intensities. A separate sample metadata
#' file is required.
#'
#' @param quant_file    Path to the matrix file (TSV, CSV, or parquet).
#' @param metadata_file Path to the sample metadata file (TSV, CSV, or Excel).
#' @param config        Named list following the ProteoForge config schema.
#'   Defaults to an empty list (uses all defaults).
#' @param feature_id_col Name of the feature ID column. If `NULL`, the first
#'   column is used.
#' @param extra_meta_cols Character vector of column names in the quant file
#'   that are feature metadata (not intensities). Auto-detected if `NULL`.
#'
#' @return A `ProteoForgeData` object.
#'
#' @examples
#' qf  <- system.file("extdata", "generic_lfq_small", "protein_matrix.tsv",
#'                    package = "proteoforge")
#' mf  <- system.file("extdata", "generic_lfq_small", "sample_metadata.tsv",
#'                    package = "proteoforge")
#' pfd <- import_generic(qf, mf)
#'
#' @export
import_generic <- function(quant_file,
                           metadata_file,
                           config           = list(),
                           feature_id_col   = NULL,
                           extra_meta_cols  = NULL) {

  stopifnot(file.exists(quant_file), file.exists(metadata_file))
  logger::log_info("import_generic: reading {quant_file}")

  # ── Read quant matrix ─────────────────────────────────────────────────────
  quant_dt <- .read_delim_auto(quant_file)

  # Identify feature ID column
  fid_col <- feature_id_col %||% names(quant_dt)[1]
  if (!fid_col %in% names(quant_dt)) {
    stop("Feature ID column '", fid_col, "' not found in quant file.", call. = FALSE)
  }

  # ── Read sample metadata ──────────────────────────────────────────────────
  smd <- .read_metadata(metadata_file)
  assert_columns(smd, c("sample_id", "condition"), "sample metadata")
  data.table::setDT(smd)

  # ── Identify sample columns (numeric columns matching metadata) ────────────
  non_sample_cols <- .infer_non_sample_cols(quant_dt, smd$sample_id,
                                            fid_col, extra_meta_cols)
  sample_cols <- setdiff(names(quant_dt), non_sample_cols)

  # Align: only keep samples present in both metadata and quant file
  shared <- intersect(smd$sample_id, sample_cols)
  missing_in_quant <- setdiff(smd$sample_id, sample_cols)
  missing_in_meta  <- setdiff(sample_cols,   smd$sample_id)

  if (length(missing_in_quant) > 0) {
    collect_warning(
      sprintf("Samples in metadata but not in quant file: %s",
              paste(missing_in_quant, collapse = ", ")),
      context = "import_generic"
    )
  }
  if (length(missing_in_meta) > 0) {
    collect_warning(
      sprintf("Quant columns not in metadata (ignored): %s",
              paste(missing_in_meta, collapse = ", ")),
      context = "import_generic"
    )
  }
  if (length(shared) == 0) {
    stop("No sample IDs match between metadata and quant file. ",
         "Check sample_id column names.", call. = FALSE)
  }

  smd <- smd[sample_id %in% shared]
  quant_sub <- quant_dt[, c(non_sample_cols, shared), with = FALSE]

  # ── Build raw_matrix ──────────────────────────────────────────────────────
  mat <- as.matrix(quant_sub[, shared, with = FALSE])
  mode(mat) <- "numeric"
  rownames(mat) <- as.character(quant_sub[[fid_col]])
  colnames(mat) <- shared

  # ── Apply zero → NA ───────────────────────────────────────────────────────
  if (isTRUE(config$preprocessing$zero_to_na %||% TRUE)) {
    n_zeros <- sum(mat == 0, na.rm = TRUE)
    if (n_zeros > 0) {
      mat[mat == 0] <- NA
      collect_warning(
        sprintf("Replaced %d zeros with NA", n_zeros),
        context = "import_generic"
      )
    }
  }

  # ── Feature metadata ──────────────────────────────────────────────────────
  feat_meta_cols <- setdiff(non_sample_cols, fid_col)
  fmd <- quant_sub[, c(fid_col, feat_meta_cols), with = FALSE]
  data.table::setnames(fmd, fid_col, "feature_id")

  # Add standard columns if missing
  if (!"protein_id"    %in% names(fmd)) fmd[, protein_id    := feature_id]
  if (!"gene_symbol"   %in% names(fmd)) {
    gene_col <- intersect(c("GeneSymbol","Gene","Genes","Gene names","Gene Symbol"),
                          names(fmd))[1]
    fmd[, gene_symbol := if (!is.na(gene_col)) fmd[[gene_col]] else NA_character_]
  }

  # ── Long format ───────────────────────────────────────────────────────────
  long <- data.table::melt(
    data.table::data.table(feature_id = rownames(mat), mat),
    id.vars      = "feature_id",
    variable.name= "sample_id",
    value.name   = "intensity",
    variable.factor = FALSE
  )

  # ── Build parameters list ─────────────────────────────────────────────────
  params <- list(
    quant_file    = quant_file,
    metadata_file = metadata_file,
    feature_id_col= fid_col,
    zero_to_na    = isTRUE(config$preprocessing$zero_to_na %||% TRUE),
    n_features    = nrow(mat),
    n_samples     = ncol(mat)
  )

  logger::log_info("import_generic: {nrow(mat)} features × {ncol(mat)} samples")

  ProteoForgeData(
    raw_long        = long,
    raw_matrix      = mat,
    feature_metadata= fmd,
    sample_metadata = smd,
    analysis_level  = config$input$analysis_level %||% "protein",
    quant_type      = "lfq",
    source_software = "generic",
    source_files    = c(quant_file, metadata_file),
    parameters      = params,
    warnings        = pf_env()$warnings
  )
}

# ── Private helpers ────────────────────────────────────────────────────────────

#' Read a delimited file adaptively (TSV, CSV, or parquet)
#' @keywords internal
.read_delim_auto <- function(path) {
  ext <- tolower(fs::path_ext(path))
  if (ext == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Package 'arrow' required to read parquet files.", call. = FALSE)
    }
    return(data.table::as.data.table(arrow::read_parquet(path)))
  }
  file_size <- fs::file_size(path)
  if (ext == "csv") {
    dt <- data.table::fread(path, sep = ",", data.table = TRUE,
                            na.strings = c("", "NA", "NaN", "Filtered", "#N/A"))
  } else {
    dt <- data.table::fread(path, sep = "\t", data.table = TRUE,
                            na.strings = c("", "NA", "NaN", "Filtered", "#N/A"))
  }
  dt
}

#' Read sample metadata from TSV, CSV, or Excel
#' @keywords internal
.read_metadata <- function(path) {
  ext <- tolower(fs::path_ext(path))
  if (ext %in% c("xlsx","xls")) {
    if (!requireNamespace("openxlsx2", quietly = TRUE)) {
      stop("Package 'openxlsx2' required to read Excel metadata.", call. = FALSE)
    }
    dt <- data.table::as.data.table(openxlsx2::read_xlsx(path))
  } else if (ext == "csv") {
    dt <- data.table::fread(path, sep = ",", data.table = TRUE)
  } else {
    dt <- data.table::fread(path, sep = "\t", data.table = TRUE)
  }
  dt
}

#' Infer non-sample columns from a quant data.table
#' @keywords internal
.infer_non_sample_cols <- function(dt, sample_ids, fid_col,
                                   extra_meta_cols = NULL) {
  # Known meta column names (from all supported software)
  known_meta_patterns <- c(
    "ProteinID","Protein\\.ID","Protein IDs","Majority protein IDs",
    "GeneSymbol","Gene\\.Symbol","Gene names","Gene Symbol","Genes","Gene",
    "Description","Protein names","Protein Description","First\\.Protein\\.Description",
    "Entry Name","Accession","Mol\\. weight","Q-value","Score",
    "Reverse","Potential contaminant","Only identified by site",
    "Master","Protein FDR Confidence","Contaminant",
    "Protein Probability","Top Peptide Probability",
    "Unique Spectral Count","Razor Spectral Count","Total Spectral Count",
    "Unique Peptide Count","Razor Peptide Count","Total Peptide Count",
    "# Peptides","# Unique Peptides","Coverage",
    "Number of proteins","Peptides","Razor.*peptides","Unique peptides",
    "Protein\\.Group","Protein\\.Ids","Protein\\.Names"
  )
  meta_pattern <- paste(known_meta_patterns, collapse = "|")

  non_sample <- unique(c(
    fid_col,
    extra_meta_cols,
    names(dt)[grepl(meta_pattern, names(dt), ignore.case = TRUE)],
    names(dt)[!sapply(dt, function(x) is.numeric(x) || all(is.na(x)))]
  ))

  # Ensure actual sample IDs are NOT in non_sample
  setdiff(non_sample, sample_ids)
}

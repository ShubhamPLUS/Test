#' Auto-detect the upstream software from file headers and names
#'
#' Examines column names and file names to identify the source software.
#' Returns a string used to dispatch to the correct importer.
#'
#' @param quant_file Path to the quantification file.
#' @param hint       Optional override: one of `"auto"`, `"diann"`,
#'   `"spectronaut"`, `"fragpipe"`, `"maxquant"`, `"pd"`, `"skyline"`,
#'   `"generic"`.
#'
#' @return A single lowercase string identifying the software.
#'
#' @examples
#' f <- system.file("extdata", "diann_example", "report.pg_matrix.tsv",
#'                  package = "proteoforge")
#' detect_software(f)
#'
#' @export
detect_software <- function(quant_file, hint = "auto") {
  stopifnot(is.character(quant_file), length(quant_file) == 1L)

  if (!identical(hint, "auto") && hint != "") {
    hint <- tolower(hint)
    valid <- c("diann","spectronaut","fragpipe","maxquant","pd","skyline","generic")
    if (hint %in% valid) return(hint)
    warning("Unknown hint '", hint, "'. Falling back to auto-detection.", call. = FALSE)
  }

  basename_lc <- tolower(fs::path_file(quant_file))
  ext         <- tolower(fs::path_ext(quant_file))

  # ── Name-based heuristics ─────────────────────────────────────────────────
  if (grepl("pg_matrix|pr_matrix|unique_genes_matrix", basename_lc)) return("diann")
  if (grepl("^report\\.(tsv|parquet|csv)$", basename_lc))            return("diann")
  if (grepl("combined_protein", basename_lc))                         return("fragpipe")
  if (grepl("combined_peptide|combined_ion", basename_lc))            return("fragpipe")
  if (grepl("proteingroups\\.txt", basename_lc))                      return("maxquant")
  if (grepl("phospho.*sites|evidence\\.txt|peptides\\.txt",
            basename_lc))                                             return("maxquant")

  # ── Header-based heuristics ───────────────────────────────────────────────
  header <- tryCatch(
    .read_header(quant_file),
    error = function(e) character(0)
  )

  if (length(header) == 0) return("generic")

  if (any(c("Protein.Group","PG.MaxLFQ","Precursor.Normalised",
             "Global.Q.Value","Lib.Q.Value") %in% header)) return("diann")

  if (any(c("EG.Qvalue","PG.Qvalue","PG.Quantity","R.FileName",
             "EG.IsDecoy") %in% header)) return("spectronaut")

  if (any(c("Protein Probability","Top Peptide Probability",
             "MaxLFQ Intensity") %in% header) ||
      any(grepl("MaxLFQ Intensity$", header)))            return("fragpipe")

  if (any(grepl("^LFQ intensity ", header)) ||
      any(grepl("^iBAQ ", header)) ||
      any(c("Majority protein IDs","Potential contaminant",
             "Only identified by site") %in% header))     return("maxquant")

  if (any(c("Master","Protein FDR Confidence: Combined",
             "Contaminant") %in% header) ||
      any(grepl("^Abundance: F[0-9]+:", header)))         return("pd")

  "generic"
}

#' Read just the header row of a delimited file
#' @keywords internal
.read_header <- function(path) {
  ext <- tolower(fs::path_ext(path))
  if (ext == "parquet") {
    if (requireNamespace("arrow", quietly = TRUE)) {
      return(names(arrow::read_parquet(path, as_data_frame = FALSE)))
    }
    return(character(0))
  }
  # TSV or CSV
  sep <- if (ext == "csv") "," else "\t"
  tryCatch({
    ln <- readLines(path, n = 1L, warn = FALSE)
    strsplit(ln, sep, fixed = TRUE)[[1]]
  }, error = function(e) character(0))
}

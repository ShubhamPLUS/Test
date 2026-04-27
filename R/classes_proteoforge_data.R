#' ProteoForgeData S4 class
#'
#' Standardized internal data container returned by all importers and carried
#' through the entire analysis pipeline.
#'
#' @slot raw_long     `data.table` — one row per feature × sample.
#' @slot raw_matrix   `matrix` — feature × sample, rownames = feature_id.
#' @slot feature_metadata `data.table` — per-feature annotations.
#' @slot sample_metadata  `data.table` — per-sample experimental design.
#' @slot design_metadata  `list` — contrasts, formula, blocking variables.
#' @slot qfeatures    `ANY` — optional QFeatures Bioconductor container.
#' @slot analysis_level `character` — one of protein/peptide/precursor/phosphosite.
#' @slot quant_type   `character` — quantification method label.
#' @slot source_software `character` — upstream software that produced the data.
#' @slot source_files `character` — paths to input files.
#' @slot parameters   `list` — full parameter list used for this object.
#' @slot parameter_hash `character` — sha256 digest of parameters.
#' @slot warnings     `character` — non-fatal messages collected during import.
#' @slot validation   `list` — validation report from metadata checks.
#'
#' @exportClass ProteoForgeData
setClass(
  "ProteoForgeData",
  representation(
    raw_long         = "ANY",
    raw_matrix       = "matrix",
    feature_metadata = "ANY",
    sample_metadata  = "ANY",
    design_metadata  = "list",
    qfeatures        = "ANY",
    analysis_level   = "character",
    quant_type       = "character",
    source_software  = "character",
    source_files     = "character",
    parameters       = "list",
    parameter_hash   = "character",
    warnings         = "character",
    validation       = "list"
  ),
  prototype(
    raw_long         = data.table::data.table(),
    raw_matrix       = matrix(nrow = 0, ncol = 0),
    feature_metadata = data.table::data.table(),
    sample_metadata  = data.table::data.table(),
    design_metadata  = list(),
    qfeatures        = NULL,
    analysis_level   = "protein",
    quant_type       = "lfq",
    source_software  = "generic",
    source_files     = character(0),
    parameters       = list(),
    parameter_hash   = character(0),
    warnings         = character(0),
    validation       = list()
  )
)

#' Validity method for ProteoForgeData
setValidity("ProteoForgeData", function(object) {
  msgs <- character(0)

  valid_levels   <- c("protein", "peptide", "precursor", "phosphosite")
  valid_quant    <- c("maxlfq", "lfq", "ibaq", "intensity", "area")
  valid_software <- c("diann", "spectronaut", "fragpipe", "maxquant",
                      "pd", "skyline", "generic")

  if (!object@analysis_level %in% valid_levels) {
    msgs <- c(msgs, sprintf(
      "analysis_level must be one of: %s", paste(valid_levels, collapse = ", ")
    ))
  }
  if (length(object@quant_type) > 0 &&
      !object@quant_type %in% valid_quant) {
    msgs <- c(msgs, sprintf(
      "quant_type must be one of: %s", paste(valid_quant, collapse = ", ")
    ))
  }
  if (!object@source_software %in% valid_software) {
    msgs <- c(msgs, sprintf(
      "source_software must be one of: %s", paste(valid_software, collapse = ", ")
    ))
  }

  if (length(msgs) > 0) msgs else TRUE
})

#' Constructor for ProteoForgeData
#'
#' @param raw_long     A `data.table` with columns feature_id, sample_id, intensity.
#' @param raw_matrix   A numeric matrix (features × samples).
#' @param feature_metadata A `data.table` with at minimum a `feature_id` column.
#' @param sample_metadata  A `data.table` with at minimum `sample_id`, `condition`.
#' @param design_metadata  Named list with elements `contrasts`, `formula`, `block_column`.
#' @param qfeatures    Optional QFeatures object (or NULL).
#' @param analysis_level One of `"protein"`, `"peptide"`, `"precursor"`, `"phosphosite"`.
#' @param quant_type   One of `"maxlfq"`, `"lfq"`, `"ibaq"`, `"intensity"`, `"area"`.
#' @param source_software One of `"diann"`, `"spectronaut"`, `"fragpipe"`,
#'   `"maxquant"`, `"pd"`, `"skyline"`, `"generic"`.
#' @param source_files Character vector of input file paths.
#' @param parameters   Named list of all parameters used.
#' @param parameter_hash sha256 hash of parameters (computed automatically if `""`)
#' @param warnings     Character vector of non-fatal warnings.
#' @param validation   List from metadata validation step.
#'
#' @return A validated `ProteoForgeData` object.
#'
#' @examples
#' mat <- matrix(c(1, 2, 3, 4), nrow = 2,
#'               dimnames = list(c("P1", "P2"), c("S1", "S2")))
#' smd <- data.table::data.table(sample_id = c("S1", "S2"),
#'                               condition = c("A", "B"))
#' pfd <- ProteoForgeData(raw_matrix = mat, sample_metadata = smd)
#'
#' @export
ProteoForgeData <- function(
    raw_long         = data.table::data.table(),
    raw_matrix       = matrix(nrow = 0, ncol = 0),
    feature_metadata = data.table::data.table(),
    sample_metadata  = data.table::data.table(),
    design_metadata  = list(),
    qfeatures        = NULL,
    analysis_level   = "protein",
    quant_type       = "lfq",
    source_software  = "generic",
    source_files     = character(0),
    parameters       = list(),
    parameter_hash   = "",
    warnings         = character(0),
    validation       = list()
) {
  if (nchar(parameter_hash) == 0 && length(parameters) > 0) {
    parameter_hash <- hash_params(parameters)
  }
  methods::new(
    "ProteoForgeData",
    raw_long         = raw_long,
    raw_matrix       = raw_matrix,
    feature_metadata = feature_metadata,
    sample_metadata  = sample_metadata,
    design_metadata  = design_metadata,
    qfeatures        = qfeatures,
    analysis_level   = analysis_level,
    quant_type       = quant_type,
    source_software  = source_software,
    source_files     = source_files,
    parameters       = parameters,
    parameter_hash   = parameter_hash,
    warnings         = warnings,
    validation       = validation
  )
}

#' Show method for ProteoForgeData
#'
#' @param object A `ProteoForgeData` object.
setMethod("show", "ProteoForgeData", function(object) {
  cat("ProteoForgeData\n")
  cat("  Software:       ", object@source_software, "\n")
  cat("  Analysis level: ", object@analysis_level, "\n")
  cat("  Quant type:     ", object@quant_type, "\n")
  cat("  Features:       ", nrow(object@raw_matrix), "\n")
  cat("  Samples:        ", ncol(object@raw_matrix), "\n")
  cat("  Warnings:       ", length(object@warnings), "\n")
  cat("  Parameter hash: ", object@parameter_hash, "\n")
  invisible(object)
})

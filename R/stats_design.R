#' Build a limma design matrix from sample metadata
#'
#' Constructs a model matrix using the specified formula. Sanitises condition
#' names to be valid R column names.
#'
#' @param sample_metadata A `data.table` with at minimum `sample_id` and
#'   `condition` columns.
#' @param formula_str     A character string formula, e.g.
#'   `"~ 0 + condition"` or `"~ condition + batch"`.
#'
#' @return A numeric `matrix` (design matrix) with rownames = `sample_id`.
#'
#' @examples
#' smd <- data.table::data.table(
#'   sample_id = paste0("S", 1:6),
#'   condition = rep(c("Control","Treatment"), each = 3),
#'   batch     = rep(c("b1","b2","b1"), 2)
#' )
#' build_design_matrix(smd, "~ 0 + condition")
#'
#' @export
build_design_matrix <- function(sample_metadata, formula_str = "~ 0 + condition") {
  smd <- data.table::copy(sample_metadata)

  # Sanitise factor levels for valid column names
  for (col in names(smd)) {
    if (is.character(smd[[col]])) {
      smd[[col]] <- make.names(smd[[col]])
    }
  }

  f   <- stats::as.formula(formula_str)
  dm  <- stats::model.matrix(f, data = smd)
  rownames(dm) <- sample_metadata$sample_id

  # Strip intercept name prefix from colnames when "~ 0 +"
  if (grepl("~ *0", formula_str)) {
    colnames(dm) <- gsub("^condition", "", colnames(dm))
    colnames(dm) <- make.names(colnames(dm))
  }
  dm
}

#' Build a limma contrast matrix from a list of contrast definitions
#'
#' @param contrasts    Named list of contrast definitions. Each element has
#'   `name`, `numerator`, `denominator` fields.
#' @param design       Design matrix (from `build_design_matrix()`).
#'
#' @return A numeric contrast matrix compatible with `limma::makeContrasts`.
#'
#' @examples
#' smd <- data.table::data.table(
#'   sample_id = paste0("S", 1:6),
#'   condition = rep(c("Control","Treatment"), each = 3)
#' )
#' dm <- build_design_matrix(smd, "~ 0 + condition")
#' contrasts <- list(list(name="Tx_vs_Ctrl", numerator="Treatment",
#'                        denominator="Control"))
#' build_contrast_matrix(contrasts, dm)
#'
#' @export
build_contrast_matrix <- function(contrasts, design) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' required.", call. = FALSE)
  }

  contrast_strings <- vapply(contrasts, function(ct) {
    num <- make.names(ct$numerator)
    den <- make.names(ct$denominator)
    sprintf("%s-%s", num, den)
  }, character(1))
  names(contrast_strings) <- vapply(contrasts, `[[`, character(1), "name")

  do.call(
    limma::makeContrasts,
    c(as.list(contrast_strings), list(levels = design))
  )
}

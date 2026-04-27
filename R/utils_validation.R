#' Assert that required columns exist in a data.frame/data.table
#'
#' Raises an informative error if any column is missing.
#'
#' @param dt       A `data.frame` or `data.table`.
#' @param required Character vector of required column names.
#' @param context  Label for error messages (e.g. `"sample metadata"`).
#'
#' @return Invisible `dt`.
#'
#' @examples
#' dt <- data.table::data.table(a = 1, b = 2)
#' assert_columns(dt, c("a", "b"), "test table")
#'
#' @export
assert_columns <- function(dt, required, context = "data") {
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop(sprintf(
      "[%s] Missing required columns: %s",
      context, paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(dt)
}

#' Validate that sample IDs in metadata match matrix column names
#'
#' @param sample_ids Character vector from sample metadata `sample_id` column.
#' @param matrix_cols Character vector of intensity matrix column names.
#'
#' @return Invisible `TRUE` on success; error on failure.
#'
#' @examples
#' validate_sample_id_match(c("S1", "S2"), c("S1", "S2"))
#'
#' @export
validate_sample_id_match <- function(sample_ids, matrix_cols) {
  in_meta_not_mat <- setdiff(sample_ids, matrix_cols)
  in_mat_not_meta <- setdiff(matrix_cols, sample_ids)
  msgs <- character(0)
  if (length(in_meta_not_mat) > 0) {
    msgs <- c(msgs, sprintf(
      "Sample IDs in metadata but not in matrix: %s",
      paste(in_meta_not_mat, collapse = ", ")
    ))
  }
  if (length(in_mat_not_meta) > 0) {
    msgs <- c(msgs, sprintf(
      "Matrix columns not in metadata: %s",
      paste(in_mat_not_meta, collapse = ", ")
    ))
  }
  if (length(msgs) > 0) {
    stop(paste(msgs, collapse = "\n"), call. = FALSE)
  }
  invisible(TRUE)
}

#' Check for batch-condition confounding
#'
#' Issues a hard error when every batch contains only one condition
#' (fully confounded). Warns when some batches contain only one condition.
#'
#' @param sample_metadata A `data.table` with `condition` and `batch` columns.
#'
#' @return Invisible `TRUE`; errors or warns as needed.
#'
#' @examples
#' smd <- data.table::data.table(
#'   condition = c("A","A","B","B"),
#'   batch     = c("b1","b2","b1","b2")
#' )
#' check_batch_confounding(smd)
#'
#' @export
check_batch_confounding <- function(sample_metadata) {
  if (!"batch" %in% names(sample_metadata)) return(invisible(TRUE))

  dt <- sample_metadata[, .(n_cond = data.table::uniqueN(condition)),
                         by = "batch"]
  n_confounded <- sum(dt$n_cond == 1)

  if (n_confounded == nrow(dt)) {
    stop(
      "Batch is fully confounded with condition: every batch contains ",
      "only one condition. Cannot include batch as covariate. ",
      "Consider removing the batch column or redesigning the experiment.",
      call. = FALSE
    )
  }
  if (n_confounded > 0) {
    warning(sprintf(
      "%d of %d batches contain only one condition. ",
      "Batch correction may be unreliable.",
      n_confounded, nrow(dt)
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate replicate counts per condition
#'
#' Errors when any condition has < 2 replicates; warns at < 3.
#'
#' @param sample_metadata A `data.table` with a `condition` column.
#' @param min_reps Minimum replicates required (hard error). Default 2.
#' @param warn_reps Threshold below which a warning is issued. Default 3.
#'
#' @return Invisible list with condition-level replicate counts.
#'
#' @examples
#' smd <- data.table::data.table(
#'   condition = c("A","A","A","B","B","B")
#' )
#' validate_replicate_counts(smd)
#'
#' @export
validate_replicate_counts <- function(sample_metadata,
                                      min_reps  = 2L,
                                      warn_reps = 3L) {
  counts <- sample_metadata[, .N, by = "condition"]
  under_min  <- counts[N < min_reps]
  under_warn <- counts[N < warn_reps & N >= min_reps]

  if (nrow(under_min) > 0) {
    stop(sprintf(
      "Conditions with fewer than %d replicates: %s. ",
      "Differential analysis is not valid.",
      min_reps,
      paste(sprintf("%s (n=%d)", under_min$condition, under_min$N),
            collapse = ", ")
    ), call. = FALSE)
  }
  if (nrow(under_warn) > 0) {
    warning(sprintf(
      "Conditions with fewer than %d replicates: %s. ",
      "Statistical power will be limited.",
      warn_reps,
      paste(sprintf("%s (n=%d)", under_warn$condition, under_warn$N),
            collapse = ", ")
    ), call. = FALSE)
  }
  invisible(counts)
}

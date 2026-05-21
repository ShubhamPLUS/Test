#' Validate and annotate a ProteoForgeData object's sample metadata
#'
#' Performs the following checks and operations:
#' 1. Required columns (`sample_id`, `condition`) present.
#' 2. `sample_id` matches matrix column names exactly.
#' 3. At least 2 conditions; warns at < 3 replicates, errors at < 2.
#' 4. Batch confounding check when `batch` column present.
#' 5. Paired design check when `subject_id` column present.
#' 6. Saves a validation report as JSON and (optionally) TSV.
#'
#' @param pfd        A `ProteoForgeData` object.
#' @param output_dir Optional path to write `validation_report.json`.
#'   If `NULL`, report is stored in the object only.
#'
#' @return The `pfd` object with `@validation` slot populated.
#'
#' @examples
#' qf  <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
#'                    package="proteoforge")
#' mf  <- system.file("extdata","generic_lfq_small","sample_metadata.tsv",
#'                    package="proteoforge")
#' pfd <- import_generic(qf, mf)
#' pfd <- validate_sample_metadata(pfd)
#'
#' @export
validate_sample_metadata <- function(pfd, output_dir = NULL) {
  stopifnot(methods::is(pfd, "ProteoForgeData"))

  smd    <- pfd@sample_metadata
  mat    <- pfd@raw_matrix
  report <- list(passed = TRUE, errors = character(0),
                 warnings = character(0), checks = list())

  # ── 1. Required columns ───────────────────────────────────────────────────
  tryCatch(
    assert_columns(smd, c("sample_id","condition"), "sample metadata"),
    error = function(e) {
      report$passed <<- FALSE
      report$errors <<- c(report$errors, conditionMessage(e))
    }
  )

  if (!report$passed) return(.finalise_validation(pfd, report, output_dir))

  # ── 2. sample_id ↔ matrix columns ─────────────────────────────────────────
  tryCatch(
    validate_sample_id_match(smd$sample_id, colnames(mat)),
    error = function(e) {
      report$passed <<- FALSE
      report$errors <<- c(report$errors, conditionMessage(e))
    }
  )

  # ── 3. Condition count and replicate checks ────────────────────────────────
  n_conds <- data.table::uniqueN(smd$condition)
  if (n_conds < 2) {
    report$passed  <- FALSE
    report$errors  <- c(report$errors,
                        "Need at least 2 conditions for differential analysis.")
  }
  tryCatch(
    validate_replicate_counts(smd),
    warning = function(w) report$warnings <<- c(report$warnings, conditionMessage(w)),
    error   = function(e) {
      report$passed <<- FALSE
      report$errors <<- c(report$errors, conditionMessage(e))
    }
  )

  # ── 4. Batch confounding ───────────────────────────────────────────────────
  if ("batch" %in% names(smd)) {
    tryCatch(
      check_batch_confounding(smd),
      warning = function(w) report$warnings <<- c(report$warnings, conditionMessage(w)),
      error   = function(e) {
        report$passed <<- FALSE
        report$errors <<- c(report$errors, conditionMessage(e))
      }
    )
    report$checks$batch_present     <- TRUE
    report$checks$n_batches         <- data.table::uniqueN(smd$batch)
  }

  # ── 5. Paired design check ────────────────────────────────────────────────
  if ("subject_id" %in% names(smd)) {
    paired_check <- .check_paired_design(smd)
    report$checks$paired_design    <- paired_check$valid
    if (!paired_check$valid) {
      report$warnings <- c(report$warnings, paired_check$message)
    }
  }

  # ── 6. Summary ────────────────────────────────────────────────────────────
  report$checks$n_samples    <- nrow(smd)
  report$checks$n_conditions <- n_conds
  report$checks$conditions   <- as.list(table(smd$condition))
  report$checks$sample_ids   <- smd$sample_id

  # Propagate warnings to global env
  for (w in report$warnings) collect_warning(w, "metadata_validate")
  for (e in report$errors)   {
    pf_env()$errors <- c(pf_env()$errors, e)
    logger::log_error("[metadata_validate] {e}")
  }

  .finalise_validation(pfd, report, output_dir)
}

#' @keywords internal
.finalise_validation <- function(pfd, report, output_dir) {
  pfd@validation <- report

  if (!is.null(output_dir)) {
    fs::dir_create(output_dir)
    json_path <- fs::path(output_dir, "validation_report.json")
    jsonlite::write_json(report, json_path, pretty = TRUE, auto_unbox = TRUE)
    logger::log_info("Validation report written to {json_path}")
  }

  if (!report$passed) {
    stop("Metadata validation failed:\n",
         paste(report$errors, collapse = "\n"), call. = FALSE)
  }

  logger::log_info("Metadata validation passed.")
  pfd
}

#' @keywords internal
.check_paired_design <- function(smd) {
  conds <- unique(smd$condition)
  subjects_per_cond <- lapply(conds, function(c) {
    smd[condition == c, subject_id]
  })
  names(subjects_per_cond) <- conds

  # Every subject should appear in every condition
  all_subjects <- Reduce(union, subjects_per_cond)
  missing <- lapply(subjects_per_cond, function(s) setdiff(all_subjects, s))
  n_missing <- sum(sapply(missing, length))

  if (n_missing > 0) {
    list(
      valid   = FALSE,
      message = sprintf(
        "Paired design: %d subjects are not present in all conditions. ",
        "Block-based analysis may be unreliable.",
        n_missing
      )
    )
  } else {
    list(valid = TRUE, message = "")
  }
}

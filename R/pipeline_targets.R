#' ProteoForge targets-based pipeline helpers
#'
#' Functions that wrap the pipeline steps specifically for use with the
#' `targets` workflow manager. These are thin wrappers around the core
#' functions and provide consistent return-value semantics expected by
#' `_targets.R`.
#'
#' @name pipeline_targets
NULL

#' Make output directories within a targets pipeline
#'
#' Wraps [make_output_dirs()] but caches the result in the targets store
#' so the path is consistent across targets within one run.
#'
#' @param config Pipeline config list.
#' @return Named list of directory paths (invisible).
#' @export
pf_targets_dirs <- function(config) {
  make_output_dirs(config)
}

#' Check whether targets pipeline metadata exists
#'
#' @param store Path to the targets store. Default `"_targets"`.
#' @return Logical.
#' @export
pf_targets_exists <- function(store = "_targets") {
  file.exists(store) && length(list.files(store)) > 0
}

#' Load a named target value from the store
#'
#' Convenience wrapper around `targets::tar_read()` that avoids a hard
#' `targets` dependency at load time.
#'
#' @param name  Target name (character).
#' @param store Path to the targets store. Default `"_targets"`.
#' @return The target's value, or `NULL` if not available.
#' @export
pf_tar_read <- function(name, store = "_targets") {
  if (!requireNamespace("targets", quietly = TRUE)) {
    warning("Package 'targets' not installed.")
    return(NULL)
  }
  tryCatch(
    targets::tar_read(!!rlang::sym(name), store = store),
    error = function(e) NULL
  )
}

#' Summarise targets pipeline status
#'
#' Returns a data.table of target names and their current status
#' (built / outdated / errored).
#'
#' @param store Path to the targets store. Default `"_targets"`.
#' @return A `data.table` with columns `name`, `status`, `seconds`.
#' @export
pf_tar_status <- function(store = "_targets") {
  if (!requireNamespace("targets", quietly = TRUE)) {
    warning("Package 'targets' not installed.")
    return(data.table::data.table())
  }
  tryCatch({
    df <- targets::tar_progress(store = store)
    data.table::as.data.table(df)
  }, error = function(e) data.table::data.table())
}

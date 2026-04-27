#' Initialise the ProteoForge run logger
#'
#' Creates a `logger` appender that writes to both console and a file.
#' Called once at the start of each pipeline run.
#'
#' @param log_dir Directory where `run.log` will be written.
#' @param level   Log level string: `"TRACE"`, `"DEBUG"`, `"INFO"`, `"WARN"`,
#'   `"ERROR"`.
#'
#' @return Invisible `NULL`.
#'
#' @examples
#' \dontrun{
#' init_logger(tempdir())
#' }
#'
#' @export
init_logger <- function(log_dir, level = "INFO") {
  fs::dir_create(log_dir)
  log_file <- fs::path(log_dir, "run.log")

  log_level <- switch(level,
    TRACE = logger::TRACE,
    DEBUG = logger::DEBUG,
    INFO  = logger::INFO,
    WARN  = logger::WARN,
    ERROR = logger::ERROR,
    logger::INFO
  )

  logger::log_appender(logger::appender_tee(log_file), index = 1)
  logger::log_threshold(log_level, index = 1)
  logger::log_formatter(logger::formatter_glue_or_sprintf, index = 1)
  logger::log_layout(
    logger::layout_glue_generator(
      format = "[{format(time, '%Y-%m-%d %H:%M:%S')}] [{level}] {msg}"
    ),
    index = 1
  )
  invisible(NULL)
}

#' Collect and save a warning to the warnings tracker
#'
#' Appends the message to an in-memory list stored in the ProteoForge
#' environment and optionally emits a `logger::log_warn()` call.
#'
#' @param msg      Warning message string.
#' @param context  Optional short context label (e.g. `"imputation"`).
#' @param env      Environment used for state accumulation.
#'   Defaults to `pf_env()`.
#'
#' @return Invisible `NULL`.
#'
#' @examples
#' collect_warning("More than 30% values imputed", context = "imputation")
#'
#' @export
collect_warning <- function(msg, context = "", env = pf_env()) {
  logger::log_warn("[{context}] {msg}")
  env$warnings <- c(env$warnings, sprintf("[%s] %s", context, msg))
  invisible(NULL)
}

#' Write collected warnings and errors to TSV files
#'
#' @param log_dir   Path to the `00_logs` directory.
#' @param warnings  Character vector of warning strings.
#' @param errors    Character vector of error strings (may be empty).
#'
#' @return Invisible list of paths written.
#'
#' @examples
#' \dontrun{
#' write_logs(tempdir(), c("Warning 1", "Warning 2"), character(0))
#' }
#'
#' @export
write_logs <- function(log_dir, warnings = character(0), errors = character(0)) {
  warn_path  <- fs::path(log_dir, "warnings.tsv")
  error_path <- fs::path(log_dir, "errors.tsv")

  w_dt <- data.table::data.table(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    message   = if (length(warnings) > 0) warnings else NA_character_
  )
  e_dt <- data.table::data.table(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    message   = if (length(errors) > 0) errors else NA_character_
  )

  data.table::fwrite(w_dt, warn_path, sep = "\t")
  data.table::fwrite(e_dt, error_path, sep = "\t")
  invisible(list(warnings = warn_path, errors = error_path))
}

#' Write `sessionInfo()` to a text file
#'
#' @param log_dir Directory to write `session_info.txt`.
#'
#' @return Invisible path.
#'
#' @examples
#' \dontrun{
#' write_session_info(tempdir())
#' }
#'
#' @export
write_session_info <- function(log_dir) {
  path <- fs::path(log_dir, "session_info.txt")
  writeLines(capture.output(sessionInfo()), path)
  invisible(path)
}

#' Access the package-level state environment
#'
#' Returns a persistent environment used to accumulate warnings, errors, and
#' run-level state across pipeline steps.
#'
#' @return An `environment`.
#'
#' @examples
#' env <- pf_env()
#' env$warnings
#'
#' @export
pf_env <- local({
  env <- new.env(parent = emptyenv())
  env$warnings <- character(0)
  env$errors   <- character(0)
  function() env
})

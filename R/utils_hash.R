#' Compute a sha256 hash of a parameter list
#'
#' Serialises the list to JSON (sorted keys) then computes a sha256 digest.
#' Used to tag every run for reproducibility.
#'
#' @param params A named list of parameters.
#' @return A 64-character hexadecimal sha256 string.
#'
#' @examples
#' hash_params(list(seed = 42, method = "median"))
#'
#' @export
hash_params <- function(params) {
  stopifnot(is.list(params))
  json_str <- jsonlite::toJSON(params, auto_unbox = TRUE, pretty = FALSE,
                               null = "null", force = TRUE)
  digest::digest(json_str, algo = "sha256", serialize = FALSE)
}

#' Write parameter hash to a file
#'
#' @param hash      sha256 string returned by `hash_params()`.
#' @param output_dir Path to the run output directory.
#'
#' @return Invisible path to the written file.
#'
#' @examples
#' \dontrun{
#' write_param_hash("abc123", output_dir = tempdir())
#' }
#'
#' @export
write_param_hash <- function(hash, output_dir) {
  path <- fs::path(output_dir, "parameter_hash.txt")
  writeLines(hash, path)
  invisible(path)
}

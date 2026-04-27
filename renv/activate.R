## This script is used to activate the renv project library.
## It is sourced automatically from .Rprofile when R starts.
## If renv is not installed, this is a no-op.

local({
  if (!requireNamespace("renv", quietly = TRUE)) {
    return(invisible(NULL))
  }
  renv::activate()
})

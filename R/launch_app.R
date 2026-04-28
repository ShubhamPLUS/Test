#' Launch the ProteoForge Shiny application
#'
#' Opens the interactive analysis interface in a browser window.
#'
#' @param port    Port to run the app on. Default `NULL` (auto-select).
#' @param host    Host to listen on. Default `"127.0.0.1"`.
#' @param launch  Whether to open the browser automatically. Default `TRUE`.
#' @param workers Number of async worker processes (via `future`). Default 1.
#'
#' @return Invisible; the app runs until interrupted.
#'
#' @examples
#' \dontrun{
#' launch_app()
#' launch_app(port = 3838, launch = FALSE)  # for Docker/server deployment
#' }
#'
#' @export
launch_app <- function(port    = NULL,
                        host    = "127.0.0.1",
                        launch  = TRUE,
                        workers = 1L) {
  app_dir <- system.file("shiny", package = "proteoforge")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    stop("Shiny app directory not found in proteoforge package.", call. = FALSE)
  }

  if (workers > 1L && requireNamespace("future", quietly = TRUE)) {
    future::plan(future::multisession, workers = workers)
  }

  shiny::runApp(
    appDir   = app_dir,
    port     = port,
    host     = host,
    launch.browser = launch
  )
}

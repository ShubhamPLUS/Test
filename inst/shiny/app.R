## ProteoForge Shiny Application
## Launch: shiny::runApp(system.file("shiny", package = "proteoforge"))
##    or:  proteoforge::launch_app()

library(shiny)
library(bslib)
library(promises)
library(future)

# Use multisession for async analysis
future::plan(future::multisession, workers = max(1L, parallel::detectCores() - 1L))

# Load package functions
if (!requireNamespace("proteoforge", quietly = TRUE)) {
  devtools::load_all(file.path(dirname(rstudioapi::getSourceEditorContext()$path),
                                "..", ".."), quiet = TRUE)
} else {
  library(proteoforge)
}

# Load all modules
module_dir <- file.path(dirname(sys.frame(1)$ofile), "modules")
for (f in list.files(module_dir, pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}

# ── UI ────────────────────────────────────────────────────────────────────
ui <- page_navbar(
  title = tags$span(
    tags$img(src = "logo.png", height = "28px", style = "margin-right:8px;"),
    "ProteoForge"
  ),
  theme = bs_theme(
    bootswatch = "flatly",
    primary    = "#2c7bb6",
    base_font  = font_google("Inter")
  ),
  id = "nav_panel",
  collapsible = TRUE,
  nav_panel("Welcome",         icon = icon("house"),     mod_welcome_ui("welcome")),
  nav_panel("Upload",          icon = icon("upload"),    mod_upload_ui("upload")),
  nav_panel("Column Mapping",  icon = icon("table"),     mod_colmap_ui("colmap")),
  nav_panel("Design",          icon = icon("flask"),     mod_design_ui("design")),
  nav_panel("QC",              icon = icon("chart-bar"), mod_qc_ui("qc")),
  nav_panel("Preprocessing",   icon = icon("filter"),    mod_preprocess_ui("preprocess")),
  nav_panel("Statistics",      icon = icon("calculator"),mod_stats_ui("stats")),
  nav_panel("Heatmaps",        icon = icon("th"),        mod_heatmap_ui("heatmap")),
  nav_panel("Venn / UpSet",    icon = icon("circle-nodes"), mod_venn_ui("venn")),
  nav_panel("Enrichment",      icon = icon("magnifying-glass"), mod_enrichment_ui("enrichment")),
  nav_panel("Network",         icon = icon("network-wired"),    mod_network_ui("network")),
  nav_panel("Phospho",         icon = icon("atom"),      mod_phospho_ui("phospho")),
  nav_panel("Report",          icon = icon("file-lines"),mod_report_ui("report")),
  nav_panel("Logs",            icon = icon("terminal"),  mod_logs_ui("logs")),
  nav_spacer(),
  nav_item(
    tags$a(
      icon("github"), "GitHub",
      href   = "https://github.com/proteoforge/proteoforge",
      target = "_blank"
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  # Shared reactive state
  rv <- reactiveValues(
    pfd_raw       = NULL,
    pfd_processed = NULL,
    stat_results  = NULL,
    qc_metrics    = NULL,
    enrichment    = NULL,
    network       = NULL,
    phospho       = NULL,
    config        = list(),
    dirs          = NULL,
    run_log       = character(0),
    busy          = FALSE
  )

  # Module servers — pass rv as shared state
  mod_welcome_server("welcome", rv)
  mod_upload_server("upload",   rv)
  mod_colmap_server("colmap",   rv)
  mod_design_server("design",   rv)
  mod_qc_server("qc",           rv)
  mod_preprocess_server("preprocess", rv)
  mod_stats_server("stats",     rv)
  mod_heatmap_server("heatmap", rv)
  mod_venn_server("venn",       rv)
  mod_enrichment_server("enrichment", rv)
  mod_network_server("network", rv)
  mod_phospho_server("phospho", rv)
  mod_report_server("report",   rv)
  mod_logs_server("logs",       rv)
}

shinyApp(ui = ui, server = server)

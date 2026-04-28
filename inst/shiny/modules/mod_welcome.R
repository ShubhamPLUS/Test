## Welcome module

mod_welcome_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    br(),
    fluidRow(
      column(8, offset = 2,
        card(
          card_header("Welcome to ProteoForge"),
          card_body(
            p("ProteoForge is a quantitative proteomics analysis platform supporting:"),
            tags$ul(
              tags$li("LFQ and DIA data from DIA-NN, Spectronaut, FragPipe, MaxQuant, Proteome Discoverer, Skyline"),
              tags$li("Differential abundance analysis with limma-DEqMS"),
              tags$li("Functional enrichment (ORA / GSEA) and STRING network analysis"),
              tags$li("Phosphoproteomics with kinase activity inference"),
              tags$li("Reproducible reporting (Excel, Word, PPTX, HTML)")
            ),
            hr(),
            h5("Quick Start"),
            tags$ol(
              tags$li("Upload your quantification file and sample metadata on the ", tags$b("Upload"), " tab."),
              tags$li("Map columns and define your experimental design."),
              tags$li("Run QC, preprocessing, and statistics."),
              tags$li("Download reports from the ", tags$b("Report"), " tab.")
            ),
            hr(),
            fluidRow(
              value_box(
                title    = "Supported Software",
                value    = "6",
                showcase = icon("microscope"),
                theme    = "primary"
              ),
              value_box(
                title    = "Analysis Modules",
                value    = "14",
                showcase = icon("cubes"),
                theme    = "success"
              ),
              value_box(
                title    = "Export Formats",
                value    = "5",
                showcase = icon("file-export"),
                theme    = "info"
              )
            )
          )
        )
      )
    )
  )
}

mod_welcome_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) { })
}

## Report module

mod_report_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    br(),
    fluidRow(
      column(6,
        card(
          card_header("Report Configuration"),
          card_body(
            textInput(ns("project_name"), "Project name", value = "ProteoForge Analysis"),
            textInput(ns("investigator"), "Investigator", value = ""),
            textAreaInput(ns("notes"), "Analysis notes", rows = 3),
            hr(),
            checkboxInput(ns("incl_methods"),  "Include methods section",     value = TRUE),
            checkboxInput(ns("incl_qc"),       "Include QC section",          value = TRUE),
            checkboxInput(ns("incl_enrich"),   "Include enrichment section",   value = TRUE),
            checkboxInput(ns("incl_network"),  "Include network section",      value = FALSE),
            checkboxInput(ns("incl_phospho"),  "Include phospho section",      value = FALSE)
          )
        )
      ),
      column(6,
        card(
          card_header("Generate Reports"),
          card_body(
            fluidRow(
              column(6,
                downloadButton(ns("dl_excel"), "Excel (.xlsx)",
                               class = "btn-success btn-block mb-2", icon = icon("file-excel")),
                downloadButton(ns("dl_word"),  "Word (.docx)",
                               class = "btn-primary btn-block mb-2", icon = icon("file-word")),
                downloadButton(ns("dl_pptx"),  "PowerPoint (.pptx)",
                               class = "btn-warning btn-block mb-2", icon = icon("file-powerpoint")),
                downloadButton(ns("dl_ai_json"), "AI Summary JSON",
                               class = "btn-secondary btn-block mb-2", icon = icon("robot"))
              ),
              column(6,
                actionButton(ns("render_html"), "Render HTML Report",
                             class = "btn-info btn-block mb-2", icon = icon("file-code")),
                br(),
                uiOutput(ns("html_link"))
              )
            )
          )
        )
      )
    ),
    fluidRow(
      column(12,
        card(
          card_header("Report Status"),
          card_body(uiOutput(ns("report_status")))
        )
      )
    )
  )
}

mod_report_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    get_config <- reactive({
      cfg <- rv$config
      cfg$project_name <- input$project_name
      cfg
    })

    tmp_dirs <- reactive({
      td <- file.path(tempdir(), "pf_report")
      dir.create(td, showWarnings = FALSE, recursive = TRUE)
      list(
        results_dir   = td,
        reports_dir   = td,
        de_dir        = td,
        qc_dir        = td,
        plots_dir     = td,
        enrichment_dir= td,
        networks_dir  = td
      )
    })

    output$dl_excel <- downloadHandler(
      filename = function() paste0(sanitise_name(input$project_name), "_results.xlsx"),
      content  = function(file) {
        req(rv$stat_results)
        withProgress(message = "Generating Excel...", {
          tmp <- tempfile(fileext = ".xlsx")
          tryCatch(
            export_excel(rv$stat_results, rv$qc_metrics,
                          enrichment_results = rv$enrichment,
                          dirs = tmp_dirs(), config = get_config()),
            error = function(e) showNotification(conditionMessage(e), type = "error")
          )
          xlsx_f <- list.files(tmp_dirs()$results_dir, "\\.xlsx$", full.names = TRUE)[1]
          if (!is.na(xlsx_f) && file.exists(xlsx_f)) file.copy(xlsx_f, file)
        })
      }
    )

    output$dl_word <- downloadHandler(
      filename = function() paste0(sanitise_name(input$project_name), "_report.docx"),
      content  = function(file) {
        req(rv$stat_results)
        withProgress(message = "Generating Word report...", {
          tryCatch(
            export_word(rv$stat_results, rv$qc_metrics,
                         dirs = tmp_dirs(), config = get_config()),
            error = function(e) showNotification(conditionMessage(e), type = "error")
          )
          docx_f <- list.files(tmp_dirs()$results_dir, "\\.docx$", full.names = TRUE)[1]
          if (!is.na(docx_f) && file.exists(docx_f)) file.copy(docx_f, file)
        })
      }
    )

    output$dl_pptx <- downloadHandler(
      filename = function() paste0(sanitise_name(input$project_name), "_slides.pptx"),
      content  = function(file) {
        req(rv$stat_results)
        withProgress(message = "Generating PowerPoint...", {
          tryCatch(
            export_pptx(rv$stat_results, dirs = tmp_dirs(), config = get_config()),
            error = function(e) showNotification(conditionMessage(e), type = "error")
          )
          pptx_f <- list.files(tmp_dirs()$results_dir, "\\.pptx$", full.names = TRUE)[1]
          if (!is.na(pptx_f) && file.exists(pptx_f)) file.copy(pptx_f, file)
        })
      }
    )

    output$dl_ai_json <- downloadHandler(
      filename = function() "ai_summary.json",
      content  = function(file) {
        req(rv$stat_results)
        out <- tryCatch(
          write_ai_summary_json(rv$stat_results, rv$qc_metrics %||% list(),
                                 dirs = tmp_dirs(), config = get_config()),
          error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }
        )
        if (!is.null(out) && file.exists(out)) file.copy(out, file)
      }
    )

    html_report_path <- reactiveVal(NULL)

    observeEvent(input$render_html, {
      req(rv$stat_results)
      withProgress(message = "Rendering HTML report...", {
        path <- tryCatch(
          render_quarto_report(rv$stat_results, rv$qc_metrics %||% list(),
                                enrichment_results = rv$enrichment,
                                phospho_results    = rv$phospho,
                                dirs = tmp_dirs(), config = get_config()),
          error = function(e) {
            showNotification(paste("Render failed:", conditionMessage(e)), type = "error")
            NULL
          }
        )
        html_report_path(path)
      })
    })

    output$html_link <- renderUI({
      path <- html_report_path()
      if (is.null(path)) return(NULL)
      tags$a(href = path, target = "_blank", class = "btn btn-info btn-block",
             icon("external-link"), " Open Report")
    })

    output$report_status <- renderUI({
      items <- list()
      if (!is.null(rv$pfd_raw))       items <- c(items, tags$li(icon("check", class="text-success"), " Data loaded"))
      if (!is.null(rv$pfd_processed)) items <- c(items, tags$li(icon("check", class="text-success"), " Preprocessing complete"))
      if (!is.null(rv$stat_results))  items <- c(items, tags$li(icon("check", class="text-success"), " Statistics complete"))
      if (!is.null(rv$enrichment))    items <- c(items, tags$li(icon("check", class="text-success"), " Enrichment complete"))
      if (length(items) == 0) {
        tags$p(class = "text-muted", "Run the analysis first.")
      } else {
        tags$ul(items)
      }
    })
  })
}

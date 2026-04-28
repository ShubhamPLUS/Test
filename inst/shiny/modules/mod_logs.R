## Logs module

mod_logs_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    br(),
    fluidRow(
      column(12,
        card(
          card_header(
            fluidRow(
              column(6, "Run Log"),
              column(6, style = "text-align:right;",
                actionButton(ns("refresh"), "Refresh", class = "btn-sm btn-outline-secondary",
                             icon = icon("rotate")),
                actionButton(ns("clear"),   "Clear",   class = "btn-sm btn-outline-danger ms-1",
                             icon = icon("trash"))
              )
            )
          ),
          card_body(
            tags$div(
              style = "height:400px; overflow-y:auto; background:#1e1e1e; padding:12px; border-radius:4px;",
              uiOutput(ns("log_text"))
            )
          )
        )
      )
    ),
    fluidRow(
      column(6,
        card(
          card_header("Warnings"),
          card_body(uiOutput(ns("warnings_list")))
        )
      ),
      column(6,
        card(
          card_header("Session Info"),
          card_body(
            verbatimTextOutput(ns("session_info"))
          )
        )
      )
    )
  )
}

mod_logs_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$clear, {
      rv$run_log <- character(0)
    })

    log_entries <- reactive({
      input$refresh
      c(
        rv$run_log,
        # Pull warnings from pf_env
        tryCatch({
          env <- pf_env()
          if (length(env$warnings) > 0) {
            paste0("[WARN] ", env$warnings)
          }
        }, error = function(e) character(0))
      )
    })

    output$log_text <- renderUI({
      entries <- log_entries()
      if (length(entries) == 0) {
        tags$span(style = "color:#888; font-family:monospace;",
                  "No log entries yet. Run an analysis to see output here.")
      } else {
        formatted <- lapply(rev(entries), function(e) {
          col <- if (grepl("\\[WARN\\]|\\[ERROR\\]", e)) "#f0a500" else "#7ec699"
          tags$div(style = sprintf("color:%s; font-family:monospace; font-size:0.85em;", col), e)
        })
        do.call(tagList, formatted)
      }
    })

    output$warnings_list <- renderUI({
      tryCatch({
        env <- pf_env()
        if (length(env$warnings) == 0) {
          tags$p(class = "text-success", icon("check"), " No warnings")
        } else {
          tags$ul(
            lapply(env$warnings, function(w) {
              tags$li(class = "text-warning", icon("triangle-exclamation"), " ", w)
            })
          )
        }
      }, error = function(e) tags$p("pf_env not available."))
    })

    output$session_info <- renderPrint({
      cat(sprintf("R version: %s\n", R.version$version.string))
      cat(sprintf("proteoforge: %s\n",
                   tryCatch(as.character(utils::packageVersion("proteoforge")), error = function(e) "dev")))
      cat(sprintf("shiny: %s\n",
                   tryCatch(as.character(utils::packageVersion("shiny")), error = function(e) "?")))
      cat(sprintf("data.table: %s\n",
                   tryCatch(as.character(utils::packageVersion("data.table")), error = function(e) "?")))
    })
  })
}

## Enrichment module

mod_enrichment_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    br(),
    fluidRow(
      column(3,
        card(
          card_header("Enrichment Settings"),
          card_body(
            selectInput(ns("contrast"), "Contrast", choices = NULL),
            checkboxGroupInput(ns("databases"), "Databases",
                               choices  = c("GO_BP","GO_MF","GO_CC","KEGG","Reactome"),
                               selected = c("GO_BP","KEGG","Reactome")),
            selectInput(ns("organism"), "Organism",
                        choices = c("Human" = "human", "Mouse" = "mouse",
                                    "Rat" = "rat", "Yeast" = "yeast")),
            numericInput(ns("fdr_enrich"), "FDR cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
            actionButton(ns("run_ora"),  "Run ORA",  class = "btn-primary me-2", icon = icon("play")),
            actionButton(ns("run_gsea"), "Run GSEA", class = "btn-outline-primary", icon = icon("play"))
          )
        )
      ),
      column(9,
        navset_tab(
          nav_panel("ORA",
            card(card_body(
              uiOutput(ns("ora_db_selector")),
              plotOutput(ns("ora_plot"), height = "480px"),
              DT::DTOutput(ns("ora_table"))
            ))
          ),
          nav_panel("GSEA",
            card(card_body(
              uiOutput(ns("gsea_db_selector")),
              plotOutput(ns("gsea_plot"), height = "480px"),
              DT::DTOutput(ns("gsea_table"))
            ))
          )
        )
      )
    )
  )
}

mod_enrichment_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observe({
      req(rv$stat_results)
      updateSelectInput(session, "contrast", choices = names(rv$stat_results))
    })

    ora_results  <- reactiveVal(NULL)
    gsea_results <- reactiveVal(NULL)

    observeEvent(input$run_ora, {
      req(rv$stat_results, rv$pfd_processed, input$contrast)
      rv$config$enrichment <- list(
        databases = input$databases,
        organism  = input$organism,
        fdr       = input$fdr_enrich
      )
      dt <- rv$stat_results[[input$contrast]]
      withProgress(message = "Running ORA...", {
        res <- tryCatch(
          run_ora(dt, rv$pfd_processed, rv$config),
          error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }
        )
        ora_results(res)
      })
      if (!is.null(ora_results())) showNotification("ORA complete.", type = "message")
    })

    observeEvent(input$run_gsea, {
      req(rv$stat_results, input$contrast)
      rv$config$enrichment <- list(
        databases = input$databases,
        organism  = input$organism,
        fdr       = input$fdr_enrich
      )
      dt <- rv$stat_results[[input$contrast]]
      withProgress(message = "Running GSEA...", {
        res <- tryCatch(
          run_gsea(dt, rv$config),
          error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }
        )
        gsea_results(res)
      })
      if (!is.null(gsea_results())) showNotification("GSEA complete.", type = "message")
    })

    output$ora_db_selector <- renderUI({
      req(ora_results())
      dbs <- names(ora_results())
      selectInput(ns("ora_db"), "Database", choices = dbs)
    })

    output$ora_plot <- renderPlot({
      req(ora_results(), input$ora_db)
      res <- ora_results()[[input$ora_db]]
      if (is.null(res) || nrow(res) == 0) { plot.new(); text(0.5, 0.5, "No results"); return() }
      tryCatch(plot_ora_dotplot(res, title = input$ora_db), error = function(e) NULL)
    })

    output$ora_table <- DT::renderDT({
      req(ora_results(), input$ora_db)
      res <- ora_results()[[input$ora_db]]
      if (is.null(res)) return(DT::datatable(data.frame()))
      DT::datatable(res, options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)
    })

    output$gsea_db_selector <- renderUI({
      req(gsea_results())
      selectInput(ns("gsea_db"), "Database", choices = names(gsea_results()))
    })

    output$gsea_plot <- renderPlot({
      req(gsea_results(), input$gsea_db)
      res <- gsea_results()[[input$gsea_db]]
      if (is.null(res) || nrow(res) == 0) { plot.new(); text(0.5, 0.5, "No results"); return() }
      tryCatch(plot_gsea_ridge(res, title = input$gsea_db), error = function(e) NULL)
    })

    output$gsea_table <- DT::renderDT({
      req(gsea_results(), input$gsea_db)
      res <- gsea_results()[[input$gsea_db]]
      if (is.null(res)) return(DT::datatable(data.frame()))
      DT::datatable(res, options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)
    })
  })
}

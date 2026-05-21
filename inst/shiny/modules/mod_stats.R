## Statistics module

mod_stats_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    br(),
    fluidRow(
      column(4,
        card(
          card_header("DE Parameters"),
          card_body(
            numericInput(ns("fdr"),   "FDR threshold",     value = 0.05, min = 0, max = 1, step = 0.01),
            numericInput(ns("log2fc"),"log₂FC threshold",  value = 1.0,  min = 0, step = 0.1),
            checkboxInput(ns("use_deqms"), "Use DEqMS (peptide-count moderation)", value = TRUE),
            actionButton(ns("run_stats"), "Run Statistics",
                         class = "btn-primary", icon = icon("play"))
          )
        )
      ),
      column(8,
        card(
          card_header("Results Summary"),
          card_body(
            uiOutput(ns("stats_badges")),
            br(),
            DT::DTOutput(ns("summary_table"))
          )
        )
      )
    ),
    fluidRow(
      column(12,
        card(
          card_header("Volcano Plot"),
          card_body(
            uiOutput(ns("contrast_selector")),
            plotOutput(ns("volcano_plot"), height = "450px")
          )
        )
      )
    ),
    fluidRow(
      column(12,
        card(
          card_header("Results Table"),
          card_body(
            uiOutput(ns("contrast_table_selector")),
            DT::DTOutput(ns("results_table")),
            br(),
            downloadButton(ns("dl_results"), "Download TSV", class = "btn-sm btn-outline-primary")
          )
        )
      )
    )
  )
}

mod_stats_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$run_stats, {
      req(rv$pfd_processed)
      rv$config$stats <- list(
        fdr           = input$fdr,
        log2fc_cutoff = input$log2fc,
        use_deqms     = input$use_deqms
      )
      withProgress(message = "Running limma-DEqMS...", {
        dm <- tryCatch(build_design_matrix(rv$pfd_processed, rv$config),
                       error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
        req(dm)
        cm <- tryCatch(build_contrast_matrix(dm, rv$config),
                       error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
        req(cm)
        rv$stat_results <- tryCatch(
          run_limma_deqms(rv$pfd_processed, dm, cm, rv$config),
          error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }
        )
      })
      if (!is.null(rv$stat_results)) {
        showNotification("Statistics complete.", type = "message")
      }
    })

    output$stats_badges <- renderUI({
      req(rv$stat_results)
      lapply(names(rv$stat_results), function(cname) {
        dt  <- rv$stat_results[[cname]]
        n   <- if ("significant" %in% names(dt)) sum(dt$significant, na.rm = TRUE) else 0L
        tags$span(class = "badge bg-primary me-2",
                   sprintf("%s: %d sig.", cname, n))
      })
    })

    output$summary_table <- DT::renderDT({
      req(rv$stat_results)
      rows <- lapply(names(rv$stat_results), function(cname) {
        dt <- rv$stat_results[[cname]]
        data.table::data.table(
          Contrast   = cname,
          Tested     = nrow(dt),
          Significant= if ("significant" %in% names(dt)) sum(dt$significant, na.rm = TRUE) else NA_integer_,
          Up         = if ("direction" %in% names(dt)) sum(dt$direction == "up", na.rm = TRUE) else NA_integer_,
          Down       = if ("direction" %in% names(dt)) sum(dt$direction == "down", na.rm = TRUE) else NA_integer_
        )
      })
      DT::datatable(data.table::rbindlist(rows), rownames = FALSE,
                    options = list(dom = "t"))
    })

    output$contrast_selector <- renderUI({
      req(rv$stat_results)
      selectInput(ns("selected_contrast"), "Select contrast",
                  choices = names(rv$stat_results), width = "300px")
    })

    output$volcano_plot <- renderPlot({
      req(rv$stat_results, input$selected_contrast)
      dt <- rv$stat_results[[input$selected_contrast]]
      tryCatch(
        plot_volcano(dt, contrast = input$selected_contrast),
        error = function(e) { plot.new(); text(0.5, 0.5, conditionMessage(e)) }
      )
    })

    output$contrast_table_selector <- renderUI({
      req(rv$stat_results)
      selectInput(ns("tbl_contrast"), "Contrast",
                  choices = names(rv$stat_results), width = "300px")
    })

    output$results_table <- DT::renderDT({
      req(rv$stat_results, input$tbl_contrast)
      dt <- rv$stat_results[[input$tbl_contrast]]
      DT::datatable(dt, options = list(scrollX = TRUE, pageLength = 15), rownames = FALSE) |>
        DT::formatRound(columns = intersect(c("log2FC","adj.P.Val","t"), names(dt)), digits = 4)
    })

    output$dl_results <- downloadHandler(
      filename = function() paste0(input$tbl_contrast, "_results.tsv"),
      content  = function(file) {
        req(rv$stat_results, input$tbl_contrast)
        data.table::fwrite(rv$stat_results[[input$tbl_contrast]], file, sep = "\t")
      }
    )
  })
}

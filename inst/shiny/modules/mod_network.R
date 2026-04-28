## Network module

mod_network_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    br(),
    fluidRow(
      column(3,
        card(
          card_header("Network Settings"),
          card_body(
            selectInput(ns("contrast"), "Contrast", choices = NULL),
            selectInput(ns("direction"), "Proteins to include",
                        choices = c("Significant (all)" = "all",
                                    "Up-regulated" = "up",
                                    "Down-regulated" = "down")),
            numericInput(ns("min_score"), "STRING min score", value = 700L, min = 0L, max = 1000L),
            numericInput(ns("top_n"),     "Top N proteins",    value = 50L,  min = 5L,  max = 500L),
            actionButton(ns("run_string"), "Run STRING",
                         class = "btn-primary", icon = icon("play"))
          )
        )
      ),
      column(9,
        card(
          card_header("STRING Network"),
          card_body(
            plotOutput(ns("network_plot"), height = "540px"),
            fluidRow(
              column(6,
                downloadButton(ns("dl_png"),   "Download PNG",    class = "btn-sm btn-outline-primary me-2"),
                downloadButton(ns("dl_graphml"),"Download GraphML", class = "btn-sm btn-outline-secondary")
              )
            )
          )
        )
      )
    ),
    fluidRow(
      column(12,
        card(
          card_header("Hub Proteins"),
          card_body(DT::DTOutput(ns("hub_table")))
        )
      )
    )
  )
}

mod_network_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observe({
      req(rv$stat_results)
      updateSelectInput(session, "contrast", choices = names(rv$stat_results))
    })

    net_r <- eventReactive(input$run_string, {
      req(rv$stat_results, input$contrast)
      dt <- rv$stat_results[[input$contrast]]
      if (input$direction != "all" && "direction" %in% names(dt)) {
        dt <- dt[direction == input$direction]
      }
      if ("significant" %in% names(dt)) dt <- dt[significant == TRUE]
      rv$config$networks <- list(
        string_min_score = input$min_score,
        top_n            = input$top_n
      )
      dirs <- list(
        results_dir  = tempdir(),
        networks_dir = tempdir()
      )
      withProgress(message = "Running STRING network...", {
        tryCatch(
          run_string_network(dt, dirs, rv$config),
          error = function(e) {
            showNotification(paste("STRING error:", conditionMessage(e)), type = "error")
            NULL
          }
        )
      })
    })

    output$network_plot <- renderPlot({
      res <- net_r()
      if (is.null(res) || is.null(res$plot)) {
        plot.new(); text(0.5, 0.5, "Network not available")
      } else {
        print(res$plot)
      }
    })

    output$hub_table <- DT::renderDT({
      res <- net_r()
      req(!is.null(res), !is.null(res$centrality))
      DT::datatable(head(res$centrality, 30),
                    options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)
    })

    output$dl_png <- downloadHandler(
      filename = function() paste0("string_network_", input$contrast, ".png"),
      content  = function(file) {
        res <- net_r()
        req(!is.null(res), !is.null(res$plot))
        ggplot2::ggsave(file, plot = res$plot, width = 10, height = 8, dpi = 150)
      }
    )

    output$dl_graphml <- downloadHandler(
      filename = function() paste0("string_network_", input$contrast, ".graphml"),
      content  = function(file) {
        res <- net_r()
        req(!is.null(res), !is.null(res$graphml_file))
        file.copy(res$graphml_file, file)
      }
    )
  })
}

## Heatmap module

mod_heatmap_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    br(),
    fluidRow(
      column(3,
        card(
          card_header("Heatmap Options"),
          card_body(
            numericInput(ns("top_n"),   "Top N features",   value = 50L,  min = 5L, max = 500L),
            selectInput(ns("sort_by"),  "Sort by",
                        choices = c("Adjusted p-value" = "adj.P.Val",
                                    "log2FC" = "log2FC",
                                    "t statistic" = "t")),
            selectInput(ns("hm_contrast"), "Contrast", choices = NULL),
            checkboxInput(ns("show_genes"), "Show gene names", value = TRUE),
            checkboxInput(ns("cluster_cols"), "Cluster columns", value = FALSE),
            actionButton(ns("draw_hm"), "Draw Heatmap",
                         class = "btn-primary", icon = icon("th"))
          )
        )
      ),
      column(9,
        card(
          card_header("Heatmap"),
          card_body(
            plotOutput(ns("heatmap_plot"), height = "600px"),
            downloadButton(ns("dl_hm"), "Download PNG", class = "btn-sm btn-outline-primary mt-2")
          )
        )
      )
    )
  )
}

mod_heatmap_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observe({
      req(rv$stat_results)
      updateSelectInput(session, "hm_contrast",
                        choices = names(rv$stat_results))
    })

    hm_plot <- eventReactive(input$draw_hm, {
      req(rv$pfd_processed, rv$stat_results, input$hm_contrast)
      dt <- rv$stat_results[[input$hm_contrast]]
      tryCatch(
        plot_heatmap(rv$pfd_processed, dt,
                     top_n        = input$top_n,
                     sort_by      = input$sort_by,
                     show_rownames= input$show_genes,
                     cluster_cols = input$cluster_cols),
        error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }
      )
    })

    output$heatmap_plot <- renderPlot({ hm_plot() })

    output$dl_hm <- downloadHandler(
      filename = function() paste0("heatmap_", input$hm_contrast, ".png"),
      content  = function(file) {
        p <- hm_plot()
        req(p)
        ggplot2::ggsave(file, plot = p, width = 10, height = 12, dpi = 150)
      }
    )
  })
}

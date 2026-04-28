## QC module

mod_qc_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    br(),
    fluidRow(
      column(12,
        card(
          card_header("Run QC Analysis"),
          card_body(
            actionButton(ns("run_qc"), "Compute QC Metrics",
                         class = "btn-primary", icon = icon("play")),
            uiOutput(ns("qc_badges"))
          )
        )
      )
    ),
    fluidRow(
      column(4,
        card(card_header("ID Counts"), card_body(plotOutput(ns("plot_id_counts"), height = "280px")))
      ),
      column(4,
        card(card_header("Intensity Distribution"), card_body(plotOutput(ns("plot_intensity"), height = "280px")))
      ),
      column(4,
        card(card_header("CV Distribution"), card_body(plotOutput(ns("plot_cv"), height = "280px")))
      )
    ),
    fluidRow(
      column(6,
        card(card_header("PCA"), card_body(plotOutput(ns("plot_pca"), height = "340px")))
      ),
      column(6,
        card(card_header("Sample Correlation"), card_body(plotOutput(ns("plot_corr"), height = "340px")))
      )
    ),
    fluidRow(
      column(12,
        card(
          card_header("Outlier Samples"),
          card_body(
            uiOutput(ns("outlier_info")),
            DT::DTOutput(ns("qc_table"))
          )
        )
      )
    )
  )
}

mod_qc_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plots_r <- reactiveVal(NULL)

    observeEvent(input$run_qc, {
      req(rv$pfd_raw)
      withProgress(message = "Computing QC metrics...", {
        rv$qc_metrics <- tryCatch(
          compute_qc_metrics(rv$pfd_raw),
          error = function(e) {
            showNotification(conditionMessage(e), type = "error")
            NULL
          }
        )
        if (!is.null(rv$qc_metrics)) {
          showNotification("QC metrics computed.", type = "message")
        }
      })
    })

    output$qc_badges <- renderUI({
      req(rv$qc_metrics)
      qc <- rv$qc_metrics
      fluidRow(
        column(3, tags$span(class = "badge bg-primary",
                             sprintf("%d samples", qc$n_samples %||% 0))),
        column(3, tags$span(class = "badge bg-info",
                             sprintf("%d features", qc$n_features %||% 0))),
        column(3, tags$span(class = "badge bg-secondary",
                             sprintf("Median CV: %.1f%%", (qc$median_cv %||% 0) * 100))),
        column(3,
          if (length(qc$outlier_samples %||% character(0)) > 0) {
            tags$span(class = "badge bg-warning",
                       sprintf("%d outlier(s)", length(qc$outlier_samples)))
          }
        )
      )
    })

    output$plot_id_counts <- renderPlot({
      req(rv$qc_metrics, rv$qc_metrics$id_counts_plot)
      rv$qc_metrics$id_counts_plot
    })
    output$plot_intensity <- renderPlot({
      req(rv$qc_metrics, rv$qc_metrics$intensity_dist_plot)
      rv$qc_metrics$intensity_dist_plot
    })
    output$plot_cv <- renderPlot({
      req(rv$qc_metrics, rv$qc_metrics$cv_plot)
      rv$qc_metrics$cv_plot
    })
    output$plot_pca <- renderPlot({
      req(rv$qc_metrics, rv$qc_metrics$pca_plot)
      rv$qc_metrics$pca_plot
    })
    output$plot_corr <- renderPlot({
      req(rv$qc_metrics, rv$qc_metrics$correlation_plot)
      rv$qc_metrics$correlation_plot
    })

    output$outlier_info <- renderUI({
      req(rv$qc_metrics)
      ol <- rv$qc_metrics$outlier_samples %||% character(0)
      if (length(ol) == 0) {
        tags$p(class = "text-success", icon("check"), " No outlier samples detected.")
      } else {
        tags$div(
          class = "alert alert-warning",
          icon("triangle-exclamation"),
          sprintf(" Outlier samples (≥2/4 criteria): %s", paste(ol, collapse = ", "))
        )
      }
    })

    output$qc_table <- DT::renderDT({
      req(rv$qc_metrics, rv$qc_metrics$per_sample_metrics)
      DT::datatable(rv$qc_metrics$per_sample_metrics,
                    options = list(pageLength = 10, scrollX = TRUE),
                    rownames = FALSE)
    })
  })
}

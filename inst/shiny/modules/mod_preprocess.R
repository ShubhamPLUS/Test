## Preprocessing module

mod_preprocess_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    br(),
    fluidRow(
      column(4,
        card(
          card_header("Filter Features"),
          card_body(
            sliderInput(ns("min_valid"), "Min valid fraction per condition",
                        min = 0, max = 1, value = 0.5, step = 0.05),
            numericInput(ns("min_reps"), "Min replicates with valid values",
                         value = 2L, min = 1L, max = 20L),
            checkboxInput(ns("rm_contam"), "Remove contaminants / reverse hits", value = TRUE)
          )
        )
      ),
      column(4,
        card(
          card_header("Normalisation"),
          card_body(
            selectInput(ns("norm_method"), "Method",
                        choices = c("Median" = "median",
                                    "Quantile" = "quantile",
                                    "VSN" = "vsn",
                                    "Cyclic loess" = "cyclic_loess",
                                    "None" = "none")),
            plotOutput(ns("before_after_norm"), height = "180px")
          )
        )
      ),
      column(4,
        card(
          card_header("Imputation"),
          card_body(
            selectInput(ns("impute_method"), "Method",
                        choices = c("Mixed MNAR/MCAR" = "mixed",
                                    "MinProb (MNAR)" = "minprob",
                                    "kNN (MCAR)" = "knn",
                                    "QRILC" = "qrilc",
                                    "Random Forest" = "rf",
                                    "None" = "none")),
            uiOutput(ns("missingness_summary"))
          )
        )
      )
    ),
    fluidRow(
      column(12,
        card(
          card_body(
            actionButton(ns("run_preprocess"), "Run Preprocessing",
                         class = "btn-primary btn-lg", icon = icon("play")),
            uiOutput(ns("pp_status"))
          )
        )
      )
    ),
    fluidRow(
      column(6,
        card(card_header("Before normalisation"), card_body(plotOutput(ns("plot_before"), height = "260px")))
      ),
      column(6,
        card(card_header("After normalisation"), card_body(plotOutput(ns("plot_after"), height = "260px")))
      )
    )
  )
}

mod_preprocess_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$run_preprocess, {
      req(rv$pfd_raw)
      rv$config$preprocessing <- list(
        min_valid_fraction = input$min_valid,
        min_replicates     = input$min_reps,
        remove_contaminants= input$rm_contam,
        normalisation      = input$norm_method,
        imputation         = input$impute_method
      )

      withProgress(message = "Preprocessing...", {
        setProgress(0.1, "Filtering features...")
        pfd <- tryCatch(filter_features(rv$pfd_raw, rv$config), error = function(e) rv$pfd_raw)
        setProgress(0.3, "Log2 transforming...")
        pfd <- log2_transform(pfd)
        setProgress(0.5, "Normalising...")
        pfd <- normalise_matrix(pfd, method = input$norm_method)
        setProgress(0.7, "Classifying missingness...")
        rv$missingness <- classify_missingness(pfd, rv$config)
        setProgress(0.85, "Imputing...")
        pfd <- tryCatch(
          impute_missing(pfd, method = input$impute_method, config = rv$config),
          error = function(e) { showNotification(conditionMessage(e), type = "warning"); pfd }
        )
        rv$pfd_processed <- pfd
        setProgress(1, "Done")
      })
      showNotification("Preprocessing complete.", type = "message")
    })

    output$pp_status <- renderUI({
      req(rv$pfd_processed)
      n <- nrow(rv$pfd_processed@raw_matrix)
      n_raw <- nrow(rv$pfd_raw@raw_matrix)
      tags$span(class = "badge bg-success ms-3",
                sprintf("%d / %d features retained", n, n_raw))
    })

    output$missingness_summary <- renderUI({
      req(rv$missingness)
      ms <- rv$missingness
      tags$div(
        tags$small(sprintf("MNAR: %d sites", ms$n_mnar %||% 0)),
        tags$br(),
        tags$small(sprintf("MCAR: %d sites", ms$n_mcar %||% 0))
      )
    })

    output$plot_before <- renderPlot({
      req(rv$pfd_raw)
      pfd <- rv$pfd_raw
      m   <- pfd@raw_matrix
      df  <- reshape2::melt(log2(m + 1))
      ggplot2::ggplot(df, ggplot2::aes(x = Var2, y = value)) +
        ggplot2::geom_boxplot(fill = "#4393c3", alpha = 0.6) +
        ggplot2::labs(x = NULL, y = "log2 Intensity") +
        theme_proteoforge() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    })

    output$plot_after <- renderPlot({
      req(rv$pfd_processed)
      m  <- rv$pfd_processed@raw_matrix
      df <- reshape2::melt(m)
      ggplot2::ggplot(df, ggplot2::aes(x = Var2, y = value)) +
        ggplot2::geom_boxplot(fill = "#74c476", alpha = 0.6) +
        ggplot2::labs(x = NULL, y = "Normalised Intensity") +
        theme_proteoforge() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    })
  })
}

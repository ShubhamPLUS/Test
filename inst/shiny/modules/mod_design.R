## Experimental design module

mod_design_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    br(),
    fluidRow(
      column(6,
        card(
          card_header("Sample Groups"),
          card_body(
            DT::DTOutput(ns("smd_table")),
            br(),
            p(class = "text-muted small",
              "Edit condition/batch assignments directly in the table, then click Apply.")
          )
        )
      ),
      column(6,
        card(
          card_header("Contrasts"),
          card_body(
            uiOutput(ns("contrast_builder")),
            br(),
            actionButton(ns("add_contrast"), "Add Contrast", icon = icon("plus"),
                         class = "btn-sm btn-outline-primary"),
            actionButton(ns("apply_design"), "Apply Design", icon = icon("check"),
                         class = "btn-success ms-2")
          )
        )
      )
    ),
    fluidRow(
      column(12,
        card(
          card_header("Design Matrix Preview"),
          card_body(verbatimTextOutput(ns("design_preview")))
        )
      )
    )
  )
}

mod_design_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Editable sample metadata table
    output$smd_table <- DT::renderDT({
      req(rv$pfd_raw)
      DT::datatable(rv$pfd_raw@sample_metadata,
                    editable  = "cell",
                    selection = "none",
                    options   = list(pageLength = 20, scrollX = TRUE),
                    rownames  = FALSE)
    })

    # Edit cell in sample metadata
    observeEvent(input$smd_table_cell_edit, {
      info <- input$smd_table_cell_edit
      req(rv$pfd_raw)
      smd <- data.table::copy(rv$pfd_raw@sample_metadata)
      smd[info$row, (names(smd)[info$col + 1]) := info$value]
      rv$pfd_raw@sample_metadata <- smd
    })

    # Dynamic contrast builder
    contrast_rows <- reactiveVal(1L)
    observeEvent(input$add_contrast, {
      contrast_rows(contrast_rows() + 1L)
    })

    output$contrast_builder <- renderUI({
      req(rv$pfd_raw)
      conditions <- unique(rv$pfd_raw@sample_metadata$condition)
      n <- contrast_rows()
      lapply(seq_len(n), function(i) {
        fluidRow(
          column(5,
            selectInput(ns(paste0("cond_a_", i)), label = if (i == 1) "Numerator" else NULL,
                        choices = conditions, width = "100%")
          ),
          column(2, tags$div(style = "text-align:center; margin-top:25px;",
                              tags$strong("vs"))),
          column(5,
            selectInput(ns(paste0("cond_b_", i)), label = if (i == 1) "Denominator" else NULL,
                        choices = conditions, width = "100%")
          )
        )
      })
    })

    observeEvent(input$apply_design, {
      req(rv$pfd_raw)
      n <- contrast_rows()
      contrasts <- character(n)
      for (i in seq_len(n)) {
        ca <- input[[paste0("cond_a_", i)]]
        cb <- input[[paste0("cond_b_", i)]]
        if (!is.null(ca) && !is.null(cb) && ca != cb) {
          contrasts[i] <- paste0(ca, "-", cb)
        }
      }
      contrasts <- contrasts[nzchar(contrasts)]
      rv$config$stats$contrasts <- contrasts
      showNotification(sprintf("%d contrast(s) set.", length(contrasts)), type = "message")
    })

    output$design_preview <- renderPrint({
      req(rv$pfd_raw)
      tryCatch({
        dm <- build_design_matrix(rv$pfd_raw, rv$config)
        print(dm)
      }, error = function(e) cat("Design matrix not yet available.\n"))
    })
  })
}

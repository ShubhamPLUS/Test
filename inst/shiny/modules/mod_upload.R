## Upload module

mod_upload_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    br(),
    fluidRow(
      column(6,
        card(
          card_header("Quantification File"),
          card_body(
            selectInput(ns("software"), "Software",
                        choices = c("Auto-detect" = "auto",
                                    "DIA-NN"       = "diann",
                                    "Spectronaut"  = "spectronaut",
                                    "FragPipe"     = "fragpipe",
                                    "MaxQuant"     = "maxquant",
                                    "Proteome Discoverer" = "pd",
                                    "Skyline"      = "skyline",
                                    "Generic TSV/CSV" = "generic")),
            fileInput(ns("quant_file"), "Quantification file",
                      accept = c(".tsv", ".csv", ".txt", ".parquet")),
            uiOutput(ns("detect_badge"))
          )
        )
      ),
      column(6,
        card(
          card_header("Sample Metadata"),
          card_body(
            fileInput(ns("meta_file"), "Sample metadata (TSV/CSV)",
                      accept = c(".tsv", ".csv", ".txt")),
            p(class = "text-muted small",
              "Required columns: sample_id, condition. Optional: batch, group, replicate."),
            downloadButton(ns("dl_template"), "Download metadata template",
                           class = "btn-sm btn-outline-secondary")
          )
        )
      )
    ),
    fluidRow(
      column(12,
        card(
          card_header("Import"),
          card_body(
            actionButton(ns("import_btn"), "Import Data",
                         class = "btn-primary", icon = icon("play")),
            tags$span(class = "ms-3"),
            uiOutput(ns("import_status"))
          )
        )
      )
    ),
    fluidRow(
      column(12,
        card(
          card_header("Data Preview"),
          card_body(
            uiOutput(ns("preview_info")),
            DT::DTOutput(ns("preview_table"))
          )
        )
      )
    )
  )
}

mod_upload_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Auto-detect badge
    output$detect_badge <- renderUI({
      req(input$quant_file)
      sw <- tryCatch(
        detect_software(input$quant_file$datapath),
        error = function(e) "unknown"
      )
      if (input$software == "auto") {
        tags$span(class = "badge bg-info ms-2",
                  paste("Detected:", sw))
      }
    })

    # Metadata template download
    output$dl_template <- downloadHandler(
      filename = "sample_metadata_template.tsv",
      content  = function(file) {
        tmpl <- data.table::data.table(
          sample_id = c("Sample_1", "Sample_2", "Sample_3", "Sample_4"),
          condition = c("Control",  "Control",  "Treatment","Treatment"),
          batch     = c("B1",       "B2",       "B1",       "B2"),
          replicate = c(1L, 2L, 1L, 2L)
        )
        data.table::fwrite(tmpl, file, sep = "\t")
      }
    )

    # Import action
    observeEvent(input$import_btn, {
      req(input$quant_file, input$meta_file)
      rv$busy <- TRUE
      shinyjs::disable("import_btn")

      p <- promises::future_promise({
        sw <- isolate(input$software)
        qf <- isolate(input$quant_file$datapath)
        mf <- isolate(input$meta_file$datapath)
        cfg <- list(input = list(software = sw))

        if (sw == "auto") sw <- detect_software(qf)
        pfd <- switch(sw,
          diann       = import_diann(qf, mf, config = cfg),
          spectronaut = import_spectronaut(qf, mf, config = cfg),
          fragpipe    = import_fragpipe(qf, mf, config = cfg),
          maxquant    = import_maxquant(qf, mf, config = cfg),
          pd          = import_pd(qf, mf, config = cfg),
          skyline     = import_skyline(qf, mf, config = cfg),
          import_generic(qf, mf, config = cfg)
        )
        pfd
      })

      p %...>% (function(pfd) {
        rv$pfd_raw <- pfd
        rv$busy    <- FALSE
        shinyjs::enable("import_btn")
        showNotification("Data imported successfully!", type = "message")
      }) %...!% (function(err) {
        rv$busy <- FALSE
        shinyjs::enable("import_btn")
        showNotification(paste("Import failed:", conditionMessage(err)), type = "error")
      })
    })

    # Status badge
    output$import_status <- renderUI({
      if (is.null(rv$pfd_raw)) return(NULL)
      n_prot <- nrow(rv$pfd_raw@raw_matrix)
      n_samp <- ncol(rv$pfd_raw@raw_matrix)
      tags$span(class = "badge bg-success",
                sprintf("%d features × %d samples", n_prot, n_samp))
    })

    # Preview
    output$preview_info <- renderUI({
      req(rv$pfd_raw)
      pfd <- rv$pfd_raw
      tags$p(class = "text-muted small",
        sprintf("Software: %s | Level: %s | Quant: %s",
                pfd@source_software, pfd@analysis_level, pfd@quant_type))
    })

    output$preview_table <- DT::renderDT({
      req(rv$pfd_raw)
      m <- rv$pfd_raw@raw_matrix
      dt <- data.table::as.data.table(m, keep.rownames = "feature_id")
      DT::datatable(head(dt, 200), options = list(scrollX = TRUE, pageLength = 10),
                    rownames = FALSE)
    })
  })
}

## Column mapping module

mod_colmap_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    br(),
    fluidRow(
      column(12,
        card(
          card_header("Column Mapping"),
          card_body(
            p(class = "text-muted",
              "Verify or correct the automatic column assignments below."),
            uiOutput(ns("colmap_ui")),
            br(),
            actionButton(ns("apply_btn"), "Apply Mapping",
                         class = "btn-success", icon = icon("check"))
          )
        )
      )
    ),
    fluidRow(
      column(12,
        card(
          card_header("Mapped Columns Preview"),
          card_body(DT::DTOutput(ns("mapped_preview")))
        )
      )
    )
  )
}

mod_colmap_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$colmap_ui <- renderUI({
      req(rv$pfd_raw)
      pfd    <- rv$pfd_raw
      fmd    <- pfd@feature_metadata
      cols   <- names(fmd)

      tagList(
        fluidRow(
          column(4,
            selectInput(ns("gene_col"), "Gene symbol column",
                        choices  = c("(none)", cols),
                        selected = intersect(c("GeneSymbol","gene_symbol","Gene"), cols)[1] %||% "(none)")
          ),
          column(4,
            selectInput(ns("prot_col"), "Protein ID column",
                        choices  = c("(none)", cols),
                        selected = intersect(c("ProteinID","protein_id","Protein"), cols)[1] %||% "(none)")
          ),
          column(4,
            selectInput(ns("feature_id_col"), "Feature ID column",
                        choices  = c("(none)", cols),
                        selected = intersect(c("feature_id","SiteID","ProteinID"), cols)[1] %||% "(none)")
          )
        ),
        fluidRow(
          column(4,
            selectInput(ns("loc_prob_col"), "Localisation probability (phospho only)",
                        choices  = c("(none)", cols),
                        selected = intersect(c("LocalizationProb","localization_prob"), cols)[1] %||% "(none)")
          ),
          column(4,
            selectInput(ns("seq_col"), "Sequence window (phospho only)",
                        choices  = c("(none)", cols),
                        selected = intersect(c("sequence_window","flanking_seq","Sequence window"), cols)[1] %||% "(none)")
          )
        )
      )
    })

    observeEvent(input$apply_btn, {
      req(rv$pfd_raw)
      pfd <- rv$pfd_raw
      fmd <- pfd@feature_metadata
      # Store mapping in config for downstream use
      rv$config$column_map <- list(
        gene_symbol      = if (input$gene_col     != "(none)") input$gene_col     else NULL,
        protein_id       = if (input$prot_col     != "(none)") input$prot_col     else NULL,
        feature_id       = if (input$feature_id_col != "(none)") input$feature_id_col else NULL,
        loc_prob         = if (input$loc_prob_col != "(none)") input$loc_prob_col else NULL,
        sequence_window  = if (input$seq_col      != "(none)") input$seq_col      else NULL
      )
      showNotification("Column mapping applied.", type = "message")
    })

    output$mapped_preview <- DT::renderDT({
      req(rv$pfd_raw)
      DT::datatable(head(rv$pfd_raw@feature_metadata, 100),
                    options = list(scrollX = TRUE, pageLength = 5),
                    rownames = FALSE)
    })
  })
}

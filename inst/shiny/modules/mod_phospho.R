## Phosphoproteomics module

mod_phospho_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    br(),
    fluidRow(
      column(3,
        card(
          card_header("Phospho Settings"),
          card_body(
            sliderInput(ns("loc_cutoff"), "Localisation probability cutoff",
                        min = 0.5, max = 1.0, value = 0.75, step = 0.05),
            checkboxInput(ns("parent_correct"), "Apply parent-protein correction", value = TRUE),
            selectInput(ns("kinase_methods"), "Kinase activity methods",
                        choices  = c("KSEA" = "ksea", "KEA3" = "kea3", "decoupleR" = "decoupler"),
                        multiple = TRUE, selected = c("ksea","decoupler")),
            numericInput(ns("consensus_min"), "Consensus min methods", value = 2L, min = 1L, max = 3L),
            br(),
            p(class = "text-muted small",
              "Note: phospho analysis requires a phosphosite-level",
              "ProteoForgeData. Use the Upload tab to load phospho data first."),
            actionButton(ns("run_phospho"), "Run Phospho Pipeline",
                         class = "btn-primary", icon = icon("atom"))
          )
        )
      ),
      column(9,
        navset_tab(
          nav_panel("DE Results",
            card(card_body(
              DT::DTOutput(ns("phospho_de_table"))
            ))
          ),
          nav_panel("Kinase Activity",
            card(card_body(
              DT::DTOutput(ns("kinase_table")),
              br(),
              DT::DTOutput(ns("consensus_table"))
            ))
          ),
          nav_panel("Volcano",
            card(card_body(
              uiOutput(ns("phospho_contrast_sel")),
              plotOutput(ns("phospho_volcano"), height = "450px")
            ))
          ),
          nav_panel("Annotation",
            card(card_body(
              DT::DTOutput(ns("annot_table"))
            ))
          )
        )
      )
    )
  )
}

mod_phospho_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$run_phospho, {
      req(rv$pfd_raw)
      if (rv$pfd_raw@analysis_level != "phosphosite") {
        showNotification("Loaded data is not phosphosite-level. Please upload a phospho dataset.",
                          type = "warning")
        return()
      }
      rv$config$phospho <- list(
        localization_prob_cutoff = input$loc_cutoff,
        parent_correct           = input$parent_correct,
        kinase_methods           = input$kinase_methods,
        consensus_min            = input$consensus_min
      )
      withProgress(message = "Running phospho pipeline...", {
        pfd_ph <- tryCatch(
          filter_phospho_sites(rv$pfd_raw, loc_cutoff = input$loc_cutoff),
          error = function(e) rv$pfd_raw
        )
        setProgress(0.2, "Preprocessing...")
        pfd_ph <- log2_transform(pfd_ph)
        pfd_ph <- normalise_matrix(pfd_ph, method = rv$config$preprocessing$normalisation %||% "median")
        pfd_ph <- tryCatch(
          impute_missing(pfd_ph, method = "minprob"),
          error = function(e) pfd_ph
        )
        setProgress(0.5, "Running statistics...")
        dm <- tryCatch(build_design_matrix(pfd_ph, rv$config), error = function(e) NULL)
        if (!is.null(dm)) {
          cm <- tryCatch(build_contrast_matrix(dm, rv$config), error = function(e) NULL)
          if (!is.null(cm)) {
            ph_stats <- tryCatch(
              run_limma_deqms(pfd_ph, dm, cm, rv$config),
              error = function(e) NULL
            )
          }
        }
        setProgress(0.7, "Kinase activity...")
        kinase_res <- if (!is.null(ph_stats) && length(ph_stats) > 0) {
          tryCatch(
            infer_kinase_activity(ph_stats[[1]],
                                   methods = input$kinase_methods,
                                   consensus_min = input$consensus_min),
            error = function(e) NULL
          )
        } else NULL

        setProgress(0.9, "Annotating...")
        annot_dt <- if (!is.null(ph_stats) && length(ph_stats) > 0) {
          tryCatch(
            annotate_phosphosites(ph_stats[[1]], use_omnipath = TRUE),
            error = function(e) ph_stats[[1]]
          )
        } else NULL

        rv$phospho <- list(
          pfd_processed = pfd_ph,
          stat_results  = ph_stats,
          kinase        = kinase_res,
          annotated     = annot_dt,
          n_sites_tested= if (!is.null(ph_stats)) nrow(ph_stats[[1]]) else 0L
        )
        setProgress(1.0, "Done")
      })
      showNotification("Phospho pipeline complete.", type = "message")
    })

    output$phospho_de_table <- DT::renderDT({
      req(rv$phospho, rv$phospho$stat_results)
      dt <- rv$phospho$stat_results[[1]]
      DT::datatable(dt, options = list(scrollX = TRUE, pageLength = 15), rownames = FALSE) |>
        DT::formatRound(columns = intersect(c("log2FC","adj.P.Val","t"), names(dt)), digits = 4)
    })

    output$kinase_table <- DT::renderDT({
      req(rv$phospho, rv$phospho$kinase)
      ks <- rv$phospho$kinase
      res <- data.table::rbindlist(
        lapply(c("ksea","kea3","decoupler"), function(m) {
          dt <- ks[[m]]
          if (is.null(dt) || nrow(dt) == 0) return(NULL)
          dt[, method := m]
          dt
        }), fill = TRUE
      )
      DT::datatable(res, options = list(scrollX = TRUE, pageLength = 15), rownames = FALSE)
    })

    output$consensus_table <- DT::renderDT({
      req(rv$phospho, rv$phospho$kinase, rv$phospho$kinase$consensus)
      DT::datatable(rv$phospho$kinase$consensus,
                    options = list(dom = "t"), rownames = FALSE)
    })

    output$phospho_contrast_sel <- renderUI({
      req(rv$phospho, rv$phospho$stat_results)
      selectInput(ns("ph_contrast"), "Contrast",
                  choices = names(rv$phospho$stat_results))
    })

    output$phospho_volcano <- renderPlot({
      req(rv$phospho, rv$phospho$stat_results, input$ph_contrast)
      dt <- rv$phospho$stat_results[[input$ph_contrast]]
      tryCatch(
        plot_volcano(dt, contrast = input$ph_contrast),
        error = function(e) { plot.new(); text(0.5, 0.5, conditionMessage(e)) }
      )
    })

    output$annot_table <- DT::renderDT({
      req(rv$phospho, rv$phospho$annotated)
      DT::datatable(rv$phospho$annotated,
                    options = list(scrollX = TRUE, pageLength = 15), rownames = FALSE)
    })
  })
}

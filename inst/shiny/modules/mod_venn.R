## Venn / UpSet module

mod_venn_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    br(),
    fluidRow(
      column(3,
        card(
          card_header("Options"),
          card_body(
            checkboxGroupInput(ns("selected_contrasts"), "Select contrasts",
                               choices = NULL),
            selectInput(ns("direction"), "Direction",
                        choices = c("Significant (all)" = "all",
                                    "Up-regulated" = "up",
                                    "Down-regulated" = "down")),
            actionButton(ns("draw_btn"), "Draw Plot",
                         class = "btn-primary", icon = icon("circle-nodes"))
          )
        )
      ),
      column(9,
        card(
          card_header("Venn / UpSet Diagram"),
          card_body(
            plotOutput(ns("venn_plot"), height = "500px"),
            downloadButton(ns("dl_venn"), "Download PNG",
                           class = "btn-sm btn-outline-primary mt-2")
          )
        )
      )
    ),
    fluidRow(
      column(12,
        card(
          card_header("Intersection Table"),
          card_body(DT::DTOutput(ns("intersect_table")))
        )
      )
    )
  )
}

mod_venn_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observe({
      req(rv$stat_results)
      updateCheckboxGroupInput(session, "selected_contrasts",
                               choices  = names(rv$stat_results),
                               selected = names(rv$stat_results))
    })

    sets_r <- eventReactive(input$draw_btn, {
      req(rv$stat_results, input$selected_contrasts)
      lapply(
        stats::setNames(input$selected_contrasts, input$selected_contrasts),
        function(cname) {
          dt <- rv$stat_results[[cname]]
          if (!"significant" %in% names(dt)) return(character(0))
          sub <- dt[significant == TRUE]
          if (input$direction != "all" && "direction" %in% names(sub)) {
            sub <- sub[direction == input$direction]
          }
          sub$feature_id
        }
      )
    })

    venn_p <- eventReactive(input$draw_btn, {
      sets <- sets_r()
      req(length(sets) >= 2)
      tryCatch(
        plot_overlap(sets),
        error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }
      )
    })

    output$venn_plot <- renderPlot({ venn_p() })

    output$dl_venn <- downloadHandler(
      filename = "venn_upset.png",
      content  = function(file) {
        p <- venn_p()
        req(p)
        ggplot2::ggsave(file, plot = p, width = 8, height = 6, dpi = 150)
      }
    )

    output$intersect_table <- DT::renderDT({
      sets <- sets_r()
      req(length(sets) >= 2)
      tryCatch({
        tbl <- get_intersection_table(sets)
        DT::datatable(tbl, options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)
      }, error = function(e) NULL)
    })
  })
}

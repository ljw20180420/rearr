library(shiny)
library(bslib)
library(tidyverse)

# Define UI ----
ui <- page_sidebar(
    title = "Micro homology",
    sidebar = sidebar(
        helpText("This APP displays the micro homology of two sequences."),
        textInput("sequence1", label = "sequence1", value = ""),
        sliderInput("seq1range", label = "Set value range", min = 0, max = 0, value = c(0, 0), step = 1),
        textInput("sequence2", label = "sequence2", value = ""),
        sliderInput("seq2range", label = "Set value range", min = 0, max = 0, value = c(0, 0), step = 1),
        numericInput("threshold", "minimal micro homology length", value = 4, min = 1)
    ),
    plotOutput("mh_matrix_plot")
)

# Define server logic ----
server <- function(input, output, session) {
    observeEvent(input$sequence1,
        {
            updateSliderInput(session, input = "seq1range",
                max = nchar(input$sequence1),
                value = c(0, nchar(input$sequence1))
            )
        }
    )
    observeEvent(input$sequence2,
        {
            updateSliderInput(session, input = "seq2range",
                max = nchar(input$sequence2),
                value = c(0, nchar(input$sequence2))
            )
        }
    )

    mh_matrix_react <- reactive({
        req(input$sequence1)
        req(input$sequence2)
        sequence1vec <- str_split(input$sequence1, "")[[1]]
        sequence2vec <- str_split(input$sequence2, "")[[1]]
        mh_matrix <- matrix(as.integer(rep(sequence1vec, time = length(sequence2vec)) == rep(sequence2vec, each = length(sequence1vec))), nrow = length(sequence1vec))
        for (i in seq_len(nrow(mh_matrix))) {
            for (j in seq_len(ncol(mh_matrix))) {
                if (i > 1 && j > 1)
                {
                if (mh_matrix[i, j] > 0) {
                    mh_matrix[i, j] <- mh_matrix[i - 1, j - 1] + mh_matrix[i, j]
                }
                }
            }
        }
        mh_matrix
    })

    pos12 <- reactive({
        req(input$threshold)
        rc <- which(mh_matrix_react() >= input$threshold, arr.ind = T)
        rc[rep(seq_len(nrow(rc)), time = mh_matrix_react()[rc]),] - rep(unlist(sapply(mh_matrix_react()[rc], function(i) seq(i - 1, 0, by = -1))), time = 2)
    })
    output$mh_matrix_plot <- renderPlot({
        tibble(pos1 = pos12()[, 1], pos2 = pos12()[, 2]) |>
            ggplot(aes(pos1, pos2)) +
            geom_point() +
            scale_x_continuous(limits=c(input$seq1range[1],input$seq1range[2]), expand=c(0,0)) +
            scale_y_continuous(limits=c(input$seq2range[1],input$seq2range[2]), expand=c(0,0)) +
            coord_equal(ratio = 1)
    })
}

# Run the app ----
shinyApp(ui = ui, server = server)
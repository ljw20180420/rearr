library(shiny)
library(bslib)
library(tidyverse)

options(shiny.maxRequestSize = 100 * 1024^2)

# Define UI ----
ui <- page_sidebar(
    title = "Base substitution frequencies",
    sidebar = sidebar(
        helpText("This APP count the frequencies of base substitution."),
        fileInput("algfile", label = "Alignment file"),
    ),
    plotOutput("baseSubFreqPlot")
)

# Define server logic ----
server <- function(input, output) {
    algLines <- reactive({
        req(input$algfile)
        readLines(input$algfile$datapath)
    })
    counts <- reactive({
        algLines()[seq(1, length(algLines()), 3)] |> strsplit("\t") |> vapply(function(x) as.integer(x[2]), 0)
    })
    subFreq <- reactive({
        refBases <- algLines()[seq(2, length(algLines()), 3)] |> toupper() |> strsplit("") |> unlist()
        queryBases <- algLines()[seq(3, length(algLines()), 3)] |> strsplit("") |> unlist()
        refLens <- vapply(algLines()[seq(2, length(algLines()), 3)], nchar, 0, USE.NAMES = FALSE)
        levels <- paste(rep(c("-", "A", "C", "G", "T"), times = 5), rep(c("-", "A", "C", "G", "T"), each = 5), sep = ">")
        tibble(sub = paste(refBases, queryBases, sep = ">"), count = rep(counts(), times = refLens)) |> summarise(count = sum(count), .by = "sub") |> mutate(sub = factor(sub, levels = levels))
    })

    output$baseSubFreqPlot <- renderPlot({
        subFreq() |>
            ggplot(aes(sub, count)) +
            geom_col() +
            scale_y_continuous(expand = c(0,0))
    })
}

# Run the app ----
shinyApp(ui = ui, server = server)
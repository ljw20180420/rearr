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
    counts <- reactive({
        req(input$algfile)
        scan(file = pipe(sprintf("sed -n '1~3p' %s | cut -f2", input$algfile$datapath)), what = integer())
    })
    refLines <- reactive({
        req(input$algfile)
        readLines(con = pipe(sprintf("sed -nr '2~3p' %s", input$algfile$datapath)))
    })
    queryLines <- reactive({
        req(input$algfile)
        readLines(con = pipe(sprintf("sed -nr '3~3p' %s", input$algfile$datapath)))
    })
    subFreq <- reactive({
        refBases <- refLines() |> toupper() |> strsplit("") |> unlist()
        queryBases <- queryLines() |> strsplit("") |> unlist()
        refLens <- vapply(refLines(), nchar, 0, USE.NAMES = FALSE)
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
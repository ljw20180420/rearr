library(shiny)
library(bslib)

# Define UI ----
ui <- page_sidebar(
    tags$head(
        tags$style(HTML("
            div.alignPair {
                font-family: Courier,courier;
                white-space: nowrap;
            }
        "))
    ),

    title = "Align one query read",
    sidebar = sidebar(
        helpText("This APP align one query read to reference1 and reference2."),
        textInput("reference1", label = "reference1", value = ""),
        textInput("reference2", label = "reference2", value = ""),
        numericInput("cut1", "cut1", value = NA, min = 0),
        numericInput("cut2", "cut2", value = NA, min = 0),
        selectInput(
            "PAM1",
            "PAM1",
            choices = c("NGG", "CCN"),
            selected = "NGG"
        ),
        selectInput(
            "PAM2",
            "PAM2",
            choices = c("NGG", "CCN"),
            selected = "NGG"
        ),
        textInput("query", label = "query read", value = "")
    ),
    htmlOutput("alignPair", class = "alignPair")
)

# Define server logic ----
server <- function(input, output, session) {
    refFile <- tempfile()
    algFile <- tempfile()
    output$alignPair <- renderText({
        req(input$reference1)
        req(input$reference2)
        req(input$cut1)
        req(input$cut2)
        req(input$query)

        writeLines(sprintf("0\n%s\n%d\n%d\n%s\n%d\n", input$reference1, input$cut1, input$cut2, input$reference2, nchar(input$reference2)), con = refFile)
        alignPipe = pipe(sprintf(
                '../../Rearrangement/build/rearrangement 3<%s -u -3 -v -9 -s0 -6 -s1 4 -s2 2 -qv -9 | gawk -f ../../tools/correct_micro_homology.AWK -- %d %s %d %s %d | tail -n+2 >%s',
                refFile,
                input$cut1,
                input$PAM1,
                input$cut2,
                input$PAM2,
                nchar(input$reference1),
                algFile
        ))
        sprintf("%s\t1", input$query) |> writeLines(con = alignPipe)
        readLines(algFile) |> paste(collapse = "<br>")
    })
}

# Run the app ----
shinyApp(ui = ui, server = server)
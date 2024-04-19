library(shiny)
library(bslib)
library(tidyverse)

# Define UI ----
ui <- page_sidebar(
    title = "Barcode Analysis",
    sidebar = sidebar(
        helpText("This APP countersides indel properties."),
        fileInput("algFiles", label = "alignment files", multiple = TRUE),
        checkboxGroupInput(
            "selectFiles",
            "Select alignment files to use",
            choices = NULL,
            selected = NULL
        )
    ),
    plotOutput("counterside")
)

# Define server logic ----
server <- function(input, output, session) {
    algTables <- vector("list", nrow(input$algFiles))
    for (i in nrow(input$algFiles)) {
        algTables[i] <- input$algFiles$datapath[i] |> sprintf("sed -n '1~3p' %s") |> read_table(
            col_names = c("index", "count", "score", "ref_id", "upDangle", "refStart1", "queryStart1", "refEnd1", "queryEnd1", "randomInsertion", "refStart2", "queryStart2", "refEnd2", "queryEnd2", "downDangle", "cut1", "cut2")
        ) |> mutate(fileName = input$algFiles$name[i])
    }

    # output$counterside <- renderPlot({

    # })
}

# Run the app ----
shinyApp(ui = ui, server = server)
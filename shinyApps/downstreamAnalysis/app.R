library(shiny)
library(bslib)
library(tidyverse)
library(shinyWidgets)

options(shiny.maxRequestSize = 100 * 1024^2)

source("helpers/alignBrowser.R")
source("helpers/alignOne.R")
source("helpers/baseSubstitute.R")
source("helpers/positionalStatistics.R")

# Define UI ----
ui <- navbarPage(
    "downstream analysis",
    tags$head(
        tags$style(HTML("
            .alignments {
                font-family: Courier,courier;
                white-space: preserve nowrap;
            }
        "))
    ),
    sidebar = sidebar(
        fileInput("algfiles", "Alignment file", multiple = TRUE),
    ),
    tabPanel("alignments browser",
        uiOutput("browserReadRangeUI"),
        htmlOutput("alignments", class = "alignments", inline = TRUE)
    ),
    tabPanel("align one pair",
        textInput("alignOneRef1", label = "reference1", value = ""),
        textInput("alignOneRef2", label = "reference2", value = ""),
        numericInput("alignOneCut1", "cut1", value = NA, min = 0),
        numericInput("alignOneCut2", "cut2", value = NA, min = 0),
        selectInput("alignOnePAM1", "PAM1", choices = c("NGG", "CCN"), selected = "NGG"),
        selectInput("alignOnePAM2", "PAM2", choices = c("NGG", "CCN"), selected = "NGG"),
        textInput("alignOneQuery", label = "query read", value = ""),
        htmlOutput("alignPair", class = "alignments", inline = TRUE)
    ),
    tabPanel("base substitution frequencies",
        plotOutput("baseSubFreqPlot")
    ),
    tabPanel("positional statistics",
        selectInput("positionalMode", "display mode", choices = c("base", "indel", "read", "snp"), selected = "base"),
        plotOutput("posBaseRef1Plot"),
        plotOutput("posBaseRef2Plot")
    )
)

# Define server logic ----
server <- function(input, output) {
    ################################
    # sidebar
    ################################
    algLines <- reactive({
        readAlgfiles(input$algfiles)
    })
    counts <- reactive({
        algLines()[seq(1, length(algLines()), 3)] |> strsplit("\t") |> vapply(function(x) as.integer(x[2]), 0)
    })

    ################################
    # alignment browser
    ################################
    output$browserReadRangeUI <- renderUI({
        req(input$algfiles)
        alignNum <- length(algLines()) / 3
        numericRangeInput("browserReadRange", "read range to browse", value = c(1, min(1000, alignNum)), min = 1, max = alignNum, step = 1)
    })
    allMd <- reactive({
        getMarkdownFromAlign(algLines(), input$browserReadRange[1], input$browserReadRange[2])
    })
    output$alignments <- renderText({
        req(input$algfiles)
        req(input$browserReadRange)
        paste(allMd(), collapse = "")
    })

    #########################################################
    # single alignment
    #########################################################
    output$alignPair <- renderText({
        req(input$alignOneRef1)
        req(input$alignOneRef2)
        req(input$alignOneCut1)
        req(input$alignOneCut2)
        req(input$alignOneQuery)

        alignOne(input$alignOneRef1, input$alignOneRef2, input$alignOneCut1, input$alignOneCut2, input$alignOnePAM1, input$alignOnePAM2, input$alignOneQuery)
    })

    #########################################################
    # base substitution frequency
    #########################################################
    subFreq <- reactive({
        countBaseSubstitute(algLines(), counts())
    })

    output$baseSubFreqPlot <- renderPlot({
        req(input$algfiles)
        subFreq() |>
            ggplot(aes(sub, count)) +
            geom_col() +
            scale_y_continuous(expand = c(0, 0))
    })

    #################################
    # positional statistics
    #################################
    cuts1 <- reactive({
        algLines()[seq(1, length(algLines()), 3)] |> strsplit("\t") |> vapply(function(x) as.integer(x[16]), 0)
    })
    ref1Lens <- reactive({
        ref1Lens <- refAllSeqs() |> vapply(function(x) gregexpr("[acgtn]", x)[[1]][2], 0)
        attr(ref1Lens, "names") <- NULL
        return(ref1Lens)
    })
    cuts2 <- reactive({
        algLines()[seq(1, length(algLines()), 3)] |> strsplit("\t") |> vapply(function(x) as.integer(x[17]), 0) - ref1Lens()
    })
    ref2Lens <- reactive({
        nchar(refAllSeqs()) - ref1Lens()
    })
    refGapMat <- reactive({
        algLines()[seq(2, length(algLines()), 3)] |> extendToSameLength()
    })
    refAllSeqs <- reactive({
        algLines()[seq(2, length(algLines()), 3)] |> str_replace_all("-", "")
    })
    ref1Mat <- reactive({
        refAllSeqs() |> substr(1, ref1Lens()) |> extendToAlignCut(cuts1())
    })
    ref2Mat <- reactive({
        refAllSeqs() |> substr(ref1Lens() + 1, nchar(refAllSeqs())) |> extendToAlignCut(cuts2())
    })
    queryAllSeqs <- reactive({
        tmpMat <- algLines()[seq(3, length(algLines()), 3)] |> extendToSameLength() |> t() |> as.vector()
        vectorToStringVector(tmpMat[(refGapMat() |> t() |> as.vector()) != '-'], nchar(refAllSeqs()))
    })
    query1Mat <- reactive({
        queryAllSeqs() |> substr(1, ref1Lens()) |> extendToAlignCut(cuts1())
    })
    query2Mat <- reactive({
        queryAllSeqs() |> substr(ref1Lens() + 1, nchar(queryAllSeqs())) |> extendToAlignCut(cuts2())
    })
    base1Tibble <- reactive({
        getPositionalBaseTibble(query1Mat(), counts(), max(cuts1()))
    })
    base2Tibble <- reactive({
        getPositionalBaseTibble(query2Mat(), counts(), max(cuts2()))
    })
    MSD1Tibble <- reactive({
        getPositionalMSDTibble(query1Mat(), ref1Mat(), counts(), max(cuts1()))
    })
    MSD2Tibble <- reactive({
        getPositionalMSDTibble(query2Mat(), ref2Mat(), counts(), max(cuts2()))
    })
    insert1Count <- reactive({
        refList <- algLines()[seq(2, length(algLines()), 3)] |> vapply(function(x) substr(x, 1, gregexpr("[acgtn]", x)[[1]][3] - 1), "") |> strsplit("")
        calInsertionCount(refList, cuts1(), ref1Lens())
    })
    insert2Count <- reactive({
        refList <- algLines()[seq(2, length(algLines()), 3)] |> vapply(function(x) substr(x, gregexpr("[acgtn]", x)[[1]][2] + 1, nchar(x)), "") |> strsplit("")
        calInsertionCount(refList, cuts2(), ref2Lens())
    })

    output$posBaseRef1Plot <- renderPlot({
        req(input$algfiles)
        req(input$positionalMode)

        if (input$positionalMode == "base") {
            drawPositionalStatic(base1Tibble(), insert1Count())
        }
        else if (input$positionalMode == "indel") {
            drawPositionalStatic(MSD1Tibble(), insert1Count())
        } else if (input$positionalMode == "read") {
            drawPositionalReads(query1Mat(), max(cuts1()))
        } else if (input$positionalMode == "snp") {
            drawPositionalSnps(query1Mat(), ref1Mat(), max(cuts1()))
        }
    })

    output$posBaseRef2Plot <- renderPlot({
        req(input$algfiles)
        req(input$positionalMode)

        if (input$positionalMode == "base") {
            drawPositionalStatic(base2Tibble(), insert2Count())
        }
        else if (input$positionalMode == "indel") {
            drawPositionalStatic(MSD2Tibble(), insert2Count())
        } else if (input$positionalMode == "read") {
            drawPositionalReads(query2Mat(), max(cuts2()))
        } else if (input$positionalMode == "snp") {
            drawPositionalSnps(query2Mat(), ref2Mat(), max(cuts2()))
        }
    })
}

# Run the app ----
shinyApp(ui = ui, server = server)
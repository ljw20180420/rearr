library(shiny)
library(bslib)
library(tidyverse)
library(ggrepel)
library(diffloop)
library(GenomicRanges)

source("../helpers/mangoFDRPValue.R")

options(shiny.maxRequestSize = 1000 * 1024^2)

ui <- page_sidebar(
    title = "Diffloop Analysis (Pair)",
    sidebar = sidebar(
        helpText("This App use diffloop to analyze loop data which require equal length."),
        fileInput("loopFiles", label = "Loop files (bedpe or mango)", multiple = TRUE),
        numericInput("mergegap", label = "merge gap", value = 500),
        numericInput("nbins", label = "nbins", value = 10, min = 0, step = 1),

        numericInput("FDR", label = "FDR", value = 1, min = 0, max = 1, step = 0.01),
        numericInput("PValue", label = "PValue", value = 1, min = 0, max = 1, step = 0.01),

        numericInput("maxgap", label = "max gap between blackList/annotation and anchor", value = 0, min = 0, step = 100),

        fileInput("blackListFiles", label = "black list regions (bed files)", multiple = TRUE),

        sliderInput("loopWidthRange", label = "loop width range", min = 0, max = 0, value = c(0, 0), step = 1),

        sliderInput("loopCountRange", label = "loop count range", min = 1, max = 1, value = c(1, 1), step = 1),

        selectInput("outputType", "output type", choices = c("loopDistancePlot", "pcaPlot"), selected = "loopDistancePlot")
    ),
    plotOutput("anyPlot")
)

server <- function(input, output) {
    beddir <- tempdir()

    samples <- reactiveVal(character(0))

    observeEvent(input$loopFiles,
        {
            if (input$outputType == "pcaPlot") {
                for (i in seq_len(length(samples()))) {
                    removeUI(selector = sprintf('.shiny-input-container:has(#%s)', samples()[i]), immediate = TRUE)
                }
            }
            samples(input$loopFiles$name |> strsplit('.', fixed = TRUE) |> vapply(function(x) x[1], ""))
            if (input$outputType == "pcaPlot") {
                for (i in seq_len(length(samples()))) {
                    insertUI(
                        selector = "div:has(>> #outputType)",
                        where = "afterEnd",
                        ui = textInput(
                            inputId = samples()[i],
                            label = samples()[i],
                            value = paste0("group", i)
                        ),
                        immediate = TRUE
                    )
                }
            }
        }
    )

    observeEvent(input$outputType, {
        if (input$outputType == "pcaPlot") {
            for (i in seq_len(length(samples()))) {
                insertUI(
                    selector = "div:has(>> #outputType)",
                    where = "afterEnd",
                    ui = textInput(
                        inputId = samples()[i],
                        label = samples()[i],
                        value = paste0("group", i)
                    ),
                    immediate = TRUE
                )
            }
        } else {
            for (i in seq_len(length(samples()))) {
                removeUI(selector = sprintf('.shiny-input-container:has(#%s)', samples()[i]), immediate = TRUE)
            }
        }
    })

    loops <- reactive({
        req(input$loopFiles)
        req(input$mergegap)
        for (i in seq_len(nrow(input$loopFiles))) {
            file.rename(input$loopFiles$datapath[i], file.path(beddir, input$loopFiles$name[i]))
        }
        if (endsWith(input$loopFiles$name[1], ".bedpe")) {
            loops <- loopsMake(beddir, samples = samples())
            loops <- mergeAnchors(loops, input$mergegap)
        } else {
            loops <- loopsMake.mango(beddir, samples = samples(), mergegap = input$mergegap)
        }
        loops |> mangoFDRPValue(nbins = input$nbins)
    })

    maskFDRPValue <- reactive({
        loops()@rowData$mango.FDR <= input$FDR & loops()@rowData$mango.P <= input$PValue
    })

    blackListRegions <- reactive({
        blackListRegions <- GRanges(c(seqnames=NULL,ranges=NULL,strand=NULL))
        if (is.null(input$blackListFiles)) {
            return(blackListRegions)
        }
        for (i in seq_len(nrow(input$blackListFiles))) {
            blackListRegionsSingle <- read_table(input$blackListFiles$datapath[i], col_names = FALSE)[,1:3] |> setNames(c("chr", "start", "end")) |> makeGRangesFromDataFrame()
            blackListRegions <- union(blackListRegions, blackListRegionsSingle)
        }
        rmchr(blackListRegions)
    })

    maskblackListLoops <- reactive({
        maskblackListAnchors <- annotateAnchors(loops(), blackListRegions(), "black", maxgap = input$maxgap)@anchors |> mcols() |> _$black
        maskblackListAnchors[loops()@interactions[, "left"]] | maskblackListAnchors[loops()@interactions[, "right"]]
    })

    observe({
        maxLoopWidth <- max(loops()@rowData$loopWidth)
        updateSliderInput(
            input = "loopWidthRange",
            max = maxLoopWidth,
            value = c(0, maxLoopWidth)
        )
    })

    maskLoopsWidth <- reactive({
        (loops()@rowData$loopWidth >= input$loopWidthRange[1]) & (loops()@rowData$loopWidth <= input$loopWidthRange[2])
    })

    observe({
        maxLoopCount <- max(max(loops()@counts), 1)
        updateSliderInput(input = "loopCountRange",
            max = maxLoopCount,
            value = c(1, maxLoopCount)
        )
    })

    maskLoopsCount <- reactive({
        loopSelectedNum <- as.matrix(loops()@counts >= input$loopCountRange[1] & loops()@counts <= input$loopCountRange[2]) |> rowSums()
        loopSelectedNum > 0
    })

    filterLoop <- reactive({
        loops() |> subsetLoops(maskFDRPValue() & !maskblackListLoops() & maskLoopsWidth() & maskLoopsCount())
    })

    output$anyPlot <- renderPlot({
        req(input$loopFiles)
        if (input$outputType == "loopDistancePlot") {
            anyPlot <- loopDistancePlot(filterLoop())
        } else if (input$outputType == "pcaPlot") {
            groups <- rep("", 3)
            loopSamples <- filterLoop()@colData |> row.names()
            for (i in seq_len(length(samples()))) {
                req(input[[samples()[i]]])
                for (j in seq_len(length(loopSamples))) {
                    if (loopSamples[j] == samples()[i]) {
                        groups[j] <- input[[samples()[i]]]
                        break
                    }
                }
            }
            anyPlot <- filterLoop() |> updateLDGroups(groups) |>
                pcaPlot() +
                geom_text_repel(aes(label = loopSamples))
        }
        return(anyPlot)
    })
}

shinyApp(ui = ui, server = server)
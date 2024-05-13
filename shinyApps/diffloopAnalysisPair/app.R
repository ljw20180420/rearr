library(shiny)
library(bslib)
library(tidyverse)
library(ggrepel)
library(diffloop)
library(GenomicRanges)
# library(Sushi) only load export functions, attach(getNamespace("Sushi")) load all functions, including those unexported ones
attach(getNamespace("Sushi"))

source("../helpers/mangoFDRPValue.R")

options(shiny.maxRequestSize = 1000 * 1024^2)

ui <- navbarPage(
    "Diffloop Analysis (Pair)",
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
    ),
    tabPanel("loopDistance",
        plotOutput("loopDistancePlot"),
        downloadButton("loopDistanceDownload")
    ),
    tabPanel("pca",
        uiOutput("sampleGroupsPca"),
        plotOutput("pcaPlot"),
        downloadButton("pcaDownload")
    ),
    tabPanel("loopArc",
        textInput("chromosome", label = "chromosome"),
        numericInput("start", label = "start", value = NA, min = 0),
        numericInput("end", label = "end", value = NA, min = 0),
        plotOutput("loopArcPlot"),
        downloadButton("loopArcDownload")
    ),
    tabPanel("diffloop",
        uiOutput("sampleGroupsDiffloop"),
        selectInput("diffloopMethod", "diffloop method", choices = c("edgeR", "Voom"), selected = "edgeR"),
        downloadButton("downloadDiffloop")
    )
)

server <- function(input, output) {
    beddir <- tempdir()

    samples <- reactive({
        input$loopFiles$name |> str_replace(".interactions.all.mango", "")
    })
    
    output$sampleGroupsPca <- renderUI({
        req(input$loopFiles)
        sampleGroups <- vector('list', nrow(input$loopFiles))
        for (i in seq_len(nrow(input$loopFiles))) {
            sampleGroups[[i]] <- textInput(inputId = paste0(samples()[i], "Pca"), label = samples()[i], value = paste0("group", i))
        }
        tagList(sampleGroups)
    })

    output$sampleGroupsDiffloop <- renderUI({
        req(input$loopFiles)
        sampleGroups <- vector('list', nrow(input$loopFiles))
        for (i in seq_len(nrow(input$loopFiles))) {
            sampleGroups[[i]] <- textInput(inputId = paste0(samples()[i], "Diffloop"), label = samples()[i], value = paste0("group", i))
        }
        tagList(sampleGroups)
    })

    loops <- reactive({
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
        req(input$loopFiles)
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
        req(input$loopFiles)
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

    loopDistanceGG <- reactive({
        ggObj <- loopDistancePlot(filterLoop())
        return(ggObj)
    })

    output$loopDistancePlot <- renderPlot({
        req(input$loopFiles)
        loopDistanceGG()
    })

    output$loopDistanceDownload <- downloadHandler(
        filename = "loopDistance.pdf",
        content = function(file) {
            ggsave(file, loopDistanceGG())
        }
    )

    pcaGG <- reactive({
        groups <- rep("", 3)
        loopSamples <- filterLoop()@colData |> row.names()
        for (i in seq_len(length(samples()))) {
            req(input[[paste0(samples()[i], "Pca")]])
            for (j in seq_len(length(loopSamples))) {
                if (loopSamples[j] == samples()[i]) {
                    groups[j] <- input[[paste0(samples()[i], "Pca")]]
                    break
                }
            }
        }
        filterLoop() |> updateLDGroups(groups) |>
            pcaPlot() +
            geom_text_repel(aes(label = loopSamples))
    })

    output$pcaPlot <- renderPlot({
        req(input$loopFiles)
        pcaGG()
    })

    output$pcaDownload <- downloadHandler(
        filename = "pca.pdf",
        content = function(file) {
            ggsave(file, pcaGG())
        }
    )

    loopArcGG <- reactive({
        region <- GRanges(seqnames = input$chromosome, ranges = IRanges(start = input$start, end = input$end))
        ggObj <- loopPlot(filterLoop(), region)
        return(ggObj)
    })

    observe({
        output$loopArcPlot <- renderPlot(
            {
                req(input$loopFiles)
                req(input$chromosome)
                req(input$start)
                req(input$end)
                loopArcGG()
            },
            height = length(samples()) * 300
        )
    })

    output$loopArcDownload <- downloadHandler(
        filename = "loopArc.pdf",
        content = function(file) {
            pdf(file)
            loopArcGG() |> replayPlot()
            dev.off()
        }
    )

    output$downloadDiffloop <- downloadHandler(
        filename = "diffloopRowData.tsv",
        content = function(file) {
            req(input$loopFiles)
            req(input$diffloopMethod)
            groups <- rep("", 3)
            loopSamples <- filterLoop()@colData |> row.names()
            for (i in seq_len(length(samples()))) {
                req(input[[paste0(samples()[i], "Diffloop")]])
                for (j in seq_len(length(loopSamples))) {
                    if (loopSamples[j] == samples()[i]) {
                        groups[j] <- input[[paste0(samples()[i], "Diffloop")]]
                        break
                    }
                }
            }
            if (input$diffloopMethod == "edgeR") {
                filterLoop() |> updateLDGroups(groups) |> quickAssoc() |> _@rowData |> write_tsv(file)
            } else {
                filterLoop() |> updateLDGroups(groups) |> quickAssocVoom() |> _@rowData |> write_tsv(file)
            }
            
        }
    )
}

shinyApp(ui = ui, server = server)
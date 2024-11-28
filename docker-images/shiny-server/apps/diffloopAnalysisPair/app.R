library(shiny)
library(bslib)
library(tidyverse)
library(ggrepel)
library(diffloop)
library(GenomicRanges)
# library(Sushi) only load exported functions, plotBedpe used in loopPlot is not exported by Sushi
plotBedpe <- getFromNamespace("plotBedpe", "Sushi")

source("../helpers/mangoFDRPValue.R")

options(shiny.maxRequestSize = 1000 * 1024^2)

ui <- navbarPage(
    "Diffloop Analysis (Pair)",
    sidebar = sidebar(
        helpText("This App use diffloop to analyze loop data which require equal length."),
        fileInput("loopFiles", label = "Loop files (bedpe or mango)", multiple = TRUE),
        textOutput("atLeastTwoFile"),
        numericInput("mergegap", label = "merge gap", value = 500),
        numericInput("nbins", label = "nbins", value = 10, min = 0, step = 1),
        numericInput("FDR", label = "FDR", value = 1, min = 0, max = 1, step = 0.01),
        numericInput("PValue", label = "PValue", value = 1, min = 0, max = 1, step = 0.01),
        numericInput("maxgap", label = "max gap between blackList/annotation and anchor", value = 0, min = 0, step = 100),
        fileInput("blackListFiles", label = "black list regions (bed files)", multiple = TRUE),
        sliderInput("loopWidthRange", label = "loop width range", min = 0, max = 0, value = c(0, 0), step = 1),
        sliderInput("loopCountRange", label = "loop count range", min = 0, max = 0, value = c(0, 0), step = 1),
    ),
    tabPanel("loopDistance",
        uiOutput("loopDistancePlot")
    ),
    tabPanel("pca",
        uiOutput("sampleGroupsPca"),
        uiOutput("pcaPlot")
    ),
    tabPanel("loopArc",
        textInput("chromosome", label = "chromosome"),
        numericInput("start", label = "start", value = NA, min = 0),
        numericInput("end", label = "end", value = NA, min = 0),
        uiOutput("loopArcPlot")
    ),
    tabPanel("diffloop",
        uiOutput("sampleGroupsDiffloop"),
        selectInput("diffloopMethod", "diffloop method", choices = c("edgeR", "Voom"), selected = "edgeR"),
        downloadButton("downloadDiffloop")
    )
)

server <- function(input, output, session) {
    ################################
    # session start
    ################################
    file.path("www", session$token, "beddir") |> dir.create(recursive=TRUE)
    ################################
    # session end
    ################################
    onSessionEnded(function() {
        unlink(file.path("www", session$token), recursive = TRUE)
    })

    output$atLeastTwoFile <- renderText({
        req(input$loopFiles)
        if (length(input$loopFiles$datapath) == 1) {
            return("Upload at least two files")
        }
        return(NULL)
    })

    samples <- reactive({
        input$loopFiles$name |> str_replace(".interactions.all.mango", "")
    })
    
    output$sampleGroupsPca <- renderUI({
        req(input$loopFiles)
        req(length(input$loopFiles$datapath) > 1)
        sampleGroups <- vector('list', nrow(input$loopFiles))
        for (i in seq_len(nrow(input$loopFiles))) {
            sampleGroups[[i]] <- textInput(inputId = paste0(samples()[i], "Pca"), label = samples()[i], value = paste0("group", i))
        }
        tagList(sampleGroups)
    })

    output$sampleGroupsDiffloop <- renderUI({
        req(input$loopFiles)
        req(length(input$loopFiles$datapath) > 1)
        sampleGroups <- vector('list', nrow(input$loopFiles))
        for (i in seq_len(nrow(input$loopFiles))) {
            sampleGroups[[i]] <- textInput(inputId = paste0(samples()[i], "Diffloop"), label = samples()[i], value = paste0("group", i))
        }
        tagList(sampleGroups)
    })

    loops <- reactive({
        for (i in seq_len(nrow(input$loopFiles))) {
            file.copy(input$loopFiles$datapath[i], file.path("www", session$token, "beddir", input$loopFiles$name[i]))
        }
        if (endsWith(input$loopFiles$name[1], ".bedpe")) {
            loops <- loopsMake(file.path("www", session$token, "beddir"), samples = samples())
            loops <- mergeAnchors(loops, input$mergegap)
        } else {
            loops <- loopsMake.mango(file.path("www", session$token, "beddir"), samples = samples(), mergegap = input$mergegap)
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

    # Use loopWidthRange(Low\High) as a proxy of input$loopWidthRange. updateSliderInput will not update input$loopWidthRange until the client send the updated values back to the server. loopWidthRange helps to block renderUI until input$loopWidthRange got updated.
    loopWidthRange <- reactiveValues(low=0, high=0)
    observe({
        loopWidthRange$low <- input$loopWidthRange[1]
        loopWidthRange$high <- input$loopWidthRange[2]
    })
    observe({
        req(input$loopFiles)
        req(length(input$loopFiles$datapath) > 1)
        loopWidthRange$low <- 0
        loopWidthRange$high <- 0
        maxLoopWidth <- max(loops()@rowData$loopWidth)
        updateSliderInput(
            input = "loopWidthRange",
            max = maxLoopWidth,
            value = c(0, maxLoopWidth)
        )
    })

    maskLoopsWidth <- reactive({
        (loops()@rowData$loopWidth >= loopWidthRange$low) & (loops()@rowData$loopWidth <= loopWidthRange$high)
    })

    # Use loopCountRange as a proxy of input$loopCountRange. updateSliderInput will not update input$loopCountRange until the client send the updated values back to the server. loopCountRange helps to block renderUI until input$loopCountRange got updated.
    loopCountRange <- reactiveValues(low=0, high=0)
    observe({
        loopCountRange$low <- input$loopCountRange[1]
        loopCountRange$high <- input$loopCountRange[2]
    })
    observe({
        req(input$loopFiles)
        req(length(input$loopFiles$datapath) > 1)
        loopCountRange$low <- 0
        loopCountRange$high <- 0
        maxLoopCount <- max(max(loops()@counts), 1)
        updateSliderInput(input = "loopCountRange",
            max = maxLoopCount,
            value = c(1, maxLoopCount)
        )
    })

    maskLoopsCount <- reactive({
        loopSelectedNum <- as.matrix(loops()@counts >= loopCountRange$low & loops()@counts <= loopCountRange$high) |> rowSums()
        loopSelectedNum > 0
    })

    filterLoop <- reactive({
        loops() |> subsetLoops(maskFDRPValue() & !maskblackListLoops() & maskLoopsWidth() & maskLoopsCount())
    })

    loopDistanceTempFile <- tempfile(tmpdir=file.path("www", session$token), fileext=".pdf")
    output$loopDistancePlot <- renderUI({
        req(input$loopFiles)
        req(length(input$loopFiles$datapath) > 1)
        req(loopWidthRange$high > 0)
        req(loopCountRange$high > 0)
        ggObj <- loopDistancePlot(filterLoop())
        ggsave(loopDistanceTempFile, plot = ggObj)
        tags$iframe(src = sub("^www/", "", loopDistanceTempFile), height = "1200px", width = "100%")
    })

    pcaTempFile <- tempfile(tmpdir=file.path("www", session$token), fileext=".pdf")
    output$pcaPlot <- renderUI({
        req(input$loopFiles)
        req(length(input$loopFiles$datapath) > 1)
        req(loopWidthRange$high > 0)
        req(loopCountRange$high > 0)
        loopSamples <- filterLoop()@colData |> row.names()
        groups <- rep("", length(loopSamples))
        for (i in seq_len(length(samples()))) {
            req(input[[paste0(samples()[i], "Pca")]])
            for (j in seq_len(length(loopSamples))) {
                if (loopSamples[j] == samples()[i]) {
                    groups[j] <- input[[paste0(samples()[i], "Pca")]]
                    break
                }
            }
        }
        ggObj <- filterLoop() |> updateLDGroups(groups) |>
            pcaPlot() +
            geom_text_repel(aes(label = loopSamples))
        ggsave(pcaTempFile, plot = ggObj)
        tags$iframe(src = sub("^www/", "", pcaTempFile), height = "1200px", width = "100%")
    })

    loopArcTempFile <- tempfile(tmpdir=file.path("www", session$token), fileext=".pdf")
    output$loopArcPlot <- renderUI({
        req(input$loopFiles)
        req(length(input$loopFiles$datapath) > 1)
        req(loopWidthRange$high > 0)
        req(loopCountRange$high > 0)
        req(input$chromosome)
        req(input$start)
        req(input$end)
        region <- GRanges(seqnames = input$chromosome, ranges = IRanges(start = input$start, end = input$end))
        pdf(loopArcTempFile)
        loopPlot(filterLoop(), region) |> replayPlot()
        dev.off()
        tags$iframe(src = sub("^www/", "", loopArcTempFile), height = "1200px", width = "100%")
    })

    output$downloadDiffloop <- downloadHandler(
        filename = "diffloopRowData.tsv",
        content = function(file) {
            req(input$loopFiles)
            req(length(input$loopFiles$datapath) > 1)
            req(input$diffloopMethod)
            loopSamples <- filterLoop()@colData |> row.names()
            groups <- rep("", length(loopSamples))
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

################################
# app end
################################
onStop(function() {
    temps <- setdiff(
        c(
            list.files(
                path="www",
                all.files=TRUE,
                full.names=TRUE,
                include.dirs=TRUE
            ),
            c(".RData", "done")
        ),
        c("www/.gitignore", "www/.", "www/..")
    )
    unlink(temps, recursive=TRUE)
})

# Run the app ----
shinyApp(ui = ui, server = server)
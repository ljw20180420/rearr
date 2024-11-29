library(shiny)
library(bslib)
library(shinyWidgets)
library(tidyverse)
library(ggrepel)
library(diffloop)
library(GenomicRanges)
# library(Sushi) only load exported functions, plotBedpe used in loopPlot is not exported by Sushi
# plotBedpe <- getFromNamespace("plotBedpe", "Sushi") cannot be used in loopPlot (maybe environment problem)
.GlobalEnv$plotBedpe <- getFromNamespace("plotBedpe", "Sushi")
# .GlobalEnv$plotBedpe <- Sushi:::plotBedpe # this cannot be used for CRAN packages

source("../helpers/mangoFDRPValue.R")

options(shiny.maxRequestSize = 1000 * 1024^2)

# This App use diffloop to analyze loop data which require equal length.
ui <- page_sidebar(
    sidebar = sidebar(
        tooltip(
            fileInput("loopFiles", label = "loop files", multiple = TRUE),
            r"(bedpe\mango loop files, at least two files are necessary)"
        ),
        tooltip(
            numericInput("mergegap", label = "merge gap", value = 500),
            "anchors within this threshold will be merged"
        ),
        tooltip(
            numericInput("nbins", label = "nbins", value = 10, min = 0, step = 1),
            "bin numbers of binomial model used to calculate loop FDR and p-value by mangoCorrection"
        ),
        tooltip(
            numericInput("FDR", label = "FDR", value = 1, min = 0, max = 1, step = 0.01),
            "FDR threshold for loops"
        ),
        tooltip(
            numericInput("PValue", label = "PValue", value = 1, min = 0, max = 1, step = 0.01),
            "p-value threshold for loops"
        ),
        tooltip(
            numericInput("maxgap", label = "max gap", value = 0, min = 0, step = 100),
            "loop anchor within this threshold to blackList/annotation is considered interacting"
        ),
        tooltip(
            fileInput("blackListFiles", label = "black list regions", multiple = TRUE),
            "bed files specifies the black list regions excluded from analysis, multiple bed files are supported"
        ),
        tooltip(
            uiOutput("globalFilters"),
            "global filters applied to all loop files"
        )
    ),
    navbarPage(
        title=NULL,
        tabPanel(
            title=tooltip(
                "loopDistance",
                "plot proportion of counts at various distances"
            ),
            tooltip(
                uiOutput("loopDistancePlot"),
                "proportion of counts at various distances"
            )
        ),
        tabPanel(
            title=tooltip(
                "pca",
                "plot principal component analysis"
            ),
            tooltip(
                uiOutput("sampleGroupsPca"),
                "sample groups"
            ),
            tooltip(
                uiOutput("pcaPlot"),
                "principal component analysis"
            )
        ),
        tabPanel(
            title=tooltip(
                "loopArc",
                "plot loop arc diagram"
            ),
            tooltip(
                textInput("chromosome", label = "chromosome"),
                "chromosome name"
            ),
            tooltip(
                numericInput("start", label = "start", value = NA, min = 0),
                "0-based chromosome start"
            ),
            tooltip(
                numericInput("end", label = "end", value = NA, min = 0),
                "0-based chromosome end"
            ),
            tooltip(
                uiOutput("loopArcPlot"),
                "loop arc diagram"
            )
        ),
        tabPanel(
            title=tooltip(
                "diffloop",
                "call different loops"
            ),
            tooltip(
                uiOutput("sampleGroupsDiffloop"),
                "sample groups"
            ),
            tooltip(
                selectInput("diffloopMethod", "diffloop method", choices = c("edgeR", "Voom"), selected = "edgeR"),
                r"(method to call different loops (edgeR\Voom))"
            ),
            tooltip(
                downloadButton("downloadDiffloop"),
                "press to download different loop file"
            )
        )
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

    samples <- reactive({
        input$loopFiles$name |> str_replace(".interactions.all.mango", "")
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

    # Use proxy$xxxxx as a proxy of input$xxxxx. update??????Input will not update input$xxxxx until the client send the updated values back to the server. proxy$xxxxx helps to block renderUI until input$xxxxx got updated.
    proxy <- reactiveValues()
    observe({
        for (name in names(proxy)) {
            proxy[[name]] <- input[[name]]
        }
    })

    output$sampleGroupsPca <- renderUI({
        req(input$loopFiles)
        req(length(input$loopFiles$datapath) > 1)
        sampleGroups <- vector('list', nrow(input$loopFiles))
        for (i in seq_len(nrow(input$loopFiles))) {
            proxy[[paste0(samples()[i], "Pca")]] <- NULL
            sampleGroups[[i]] <- textInput(inputId = paste0(samples()[i], "Pca"), label = samples()[i], value = paste0("group", i))
        }
        tagList(sampleGroups)
    })

    output$sampleGroupsDiffloop <- renderUI({
        req(input$loopFiles)
        req(length(input$loopFiles$datapath) > 1)
        sampleGroups <- vector('list', nrow(input$loopFiles))
        for (i in seq_len(nrow(input$loopFiles))) {
            proxy[[paste0(samples()[i], "Diffloop")]] <- NULL
            sampleGroups[[i]] <- textInput(inputId = paste0(samples()[i], "Diffloop"), label = samples()[i], value = paste0("group", i))
        }
        tagList(sampleGroups)
    })

    output$globalFilters <- renderUI({
        req(input$loopFiles)
        req(length(input$loopFiles$datapath) > 1)
        proxy$loopWidthRange <- NULL
        proxy$loopCountRange <- NULL
        minLoopWidth <- min(loops()@rowData$loopWidth)
        maxLoopWidth <- max(loops()@rowData$loopWidth)
        widthFilter <- numericRangeInput("loopWidthRange", "loop width range", value = c(minLoopWidth, maxLoopWidth), min = minLoopWidth, max = maxLoopWidth, step = 1)
        maxLoopCount <- max(max(loops()@counts), 2)
        minLoopCount <- min(min(loops()@counts), maxLoopCount)
        countFilter <- numericRangeInput("loopCountRange", "loop count range", value = c(max(minLoopCount, 2), maxLoopCount), min = minLoopCount, max = maxLoopCount, step = 1)
        tagList(widthFilter, countFilter)
    })

    maskLoopsWidth <- reactive({
        (loops()@rowData$loopWidth >= proxy$loopWidthRange[1]) & (loops()@rowData$loopWidth <= proxy$loopWidthRange[2])
    })

    maskLoopsCount <- reactive({
        loopSelectedNum <- as.matrix(loops()@counts >= proxy$loopCountRange[1] & loops()@counts <= proxy$loopCountRange[2]) |> rowSums()
        loopSelectedNum > 0
    })

    filterLoop <- reactive({
        loops() |> subsetLoops(maskFDRPValue() & !maskblackListLoops() & maskLoopsWidth() & maskLoopsCount())
    })

    loopDistanceTempFile <- tempfile(tmpdir=file.path("www", session$token), fileext=".pdf")
    output$loopDistancePlot <- renderUI({
        req(input$loopFiles)
        req(length(input$loopFiles$datapath) > 1)
        req(proxy$loopWidthRange)
        req(proxy$loopCountRange)
        ggObj <- loopDistancePlot(filterLoop())
        ggsave(loopDistanceTempFile, plot = ggObj)
        tags$iframe(src = sub("^www/", "", loopDistanceTempFile), height = "1200px", width = "100%")
    })

    pcaTempFile <- tempfile(tmpdir=file.path("www", session$token), fileext=".pdf")
    output$pcaPlot <- renderUI({
        req(input$loopFiles)
        req(length(input$loopFiles$datapath) > 1)
        req(proxy$loopWidthRange)
        req(proxy$loopCountRange)
        loopSamples <- filterLoop()@colData |> row.names()
        groups <- rep("", length(loopSamples))
        for (i in seq_len(length(samples()))) {
            req(proxy[[paste0(samples()[i], "Pca")]])
            for (j in seq_len(length(loopSamples))) {
                if (loopSamples[j] == samples()[i]) {
                    groups[j] <- proxy[[paste0(samples()[i], "Pca")]]
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
        req(proxy$loopWidthRange)
        req(proxy$loopCountRange)
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
                req(proxy[[paste0(samples()[i], "Diffloop")]])
                for (j in seq_len(length(loopSamples))) {
                    if (loopSamples[j] == samples()[i]) {
                        groups[j] <- proxy[[paste0(samples()[i], "Diffloop")]]
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
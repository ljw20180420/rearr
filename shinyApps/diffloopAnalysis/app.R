library(shiny)
library(bslib)
library(tidyverse)
library(diffloop)
library(GenomicRanges)
library(HiContacts)
library(rtracklayer)
library(GenomicFeatures)
library(shinyjs)
library(shinyWidgets)

chromInfo <- getChromInfoFromUCSC("hg19")

source("../helpers/mangoFDRPValue.R")
source("../helpers/stat_ecdf_weighted.R")

options(shiny.maxRequestSize = 1000 * 1024^2)

ui <- navbarPage(
    "Diffloop Analysis",
    sidebar = sidebar(
        helpText("This App use diffloop to analyze loop data."),
        fileInput("loopFiles", "Loop files (bedpe or mango)", multiple = TRUE),
        numericInput("mergegap", "merge gap", value = 500),
        numericInput("nbins", "nbins", value = 10, min = 0, step = 1),

        numericInput("FDR", "FDR", value = 1, min = 0, max = 1, step = 0.01),
        numericInput("PValue", "PValue", value = 1, min = 0, max = 1, step = 0.01),

        numericInput("maxgap", "max gap between blackList/annotation and anchor", value = 0, min = 0, step = 100),

        fileInput("blackListFiles", "black list regions (bed files)", multiple = TRUE),

        fileInput("annotFiles", "annotation files", multiple = TRUE),

        uiOutput("globalFilters"),

        fileInput("hicFiles", "hic files"),

        numericInput("filterNum", "filter number", value = 1, min = 1, step = 1),

        uiOutput("filters")
    ),
    tabPanel("loop number",
        selectInput("numberCount", "number/count", choices = c("number", "count")),
        plotOutput("loopNumPlot"),
        downloadButton("loopNumDownload")
    ),
    tabPanel("loop width",
        numericInput("loopWidthBin", "loop width bin number", value = 20, min = 0, step = 1),
        plotOutput("loopWidthPlot"),
        downloadButton("loopWidthDownload")
    ),
    tabPanel("loop width cumulative distribution",
        plotOutput("loopWidthCdfPlot"),
        downloadButton("loopWidthCdfDownload")
    ),
    tabPanel("heat map",
        checkboxInput("balanced", "balanced", value = TRUE),
        checkboxInput("horizontal", "horizontal", value = TRUE),
        selectInput("scale", "scale", choices = c('log10', 'log2', 'linear', 'exp0.2')),
        selectInput("resolution", "resolution", choices = NULL),
        selectInput("chromosome", "chromosome", choices = NULL),
        numericInput("start", "start", value = NA, min = 0),
        numericInput("end", "end", value = NA, min = 0),
        numericInput("maxDistance", "max distance", value = 1000000, min = 100000, step = 10000),
        actionButton("renderHeat", "render heatmap"),
        uiOutput("heatMapPlot")
    ),
    tabPanel("aggregate",
        useShinyjs(),
        selectInput("normalization", "normalization", choices = c("weight", "''")),
        selectInput("aggResolution", "resolution", choices = NULL),
        numericInput("flank", "flanking bins", value = 10, min = 0, step = 1),
        fileInput("view", "view file in bed which defines which regions of the chromosomes to use"),
        checkboxInput("setVminVmax", "set vmin and vmax", value = FALSE),
        numericInput("vmin", "value for the lowest colour", value = NA),
        numericInput("vmax", "value for the highest colour", value = NA),
        actionButton("renderAggregate", "render aggregate"),
        uiOutput("aggregatePlot")
    )
)

server <- function(input, output) {
    beddir <- tempdir()

    samples <- reactive({
        input$loopFiles$name |> str_replace(".interactions.all.mango", "")
    })

    resolutions <- reactive({
        resolutions <- availableResolutions(input$hicFiles$datapath)
        return(resolutions[resolutions >= 10000])
    })

    observe({
        toggleState("vmin")
        toggleState("vmax")
    }) |> bindEvent(input$setVminVmax)

    observe({
        req(input$hicFiles)
        updateSelectInput(inputId = "resolution", choices = resolutions())
        updateSelectInput(inputId = "aggResolution", choices = resolutions())
        updateSelectInput(inputId = "chromosome", choices = availableChromosomes(input$hicFiles$datapath)@seqnames)
    })

    observe({
        req(input$chromosome)
        chromSize <- chromInfo$size[chromInfo$chrom == input$chromosome]
        updateNumericInput(inputId = "start", value = 0)
        updateNumericInput(inputId = "end", value = chromSize, max = chromSize)
    })

    output$filters <- renderUI({
        req(input$loopFiles)
        annoteNum <- ifelse(is.null(input$annotFiles), 0, nrow(input$annotFiles))
        filterUIs <- vector('list', input$filterNum * (annoteNum + 2))
        for (i in seq_len(input$filterNum)) {
            for (j in seq_len(annoteNum + 2)) {
                p <- (i - 1) * (annoteNum + 2) + j
                if (j == 1) {
                    filterUIs[[p]] <- checkboxGroupInput(paste0("filter", i), paste0("filter", i), choices = samples(), selected = samples())
                } else if (j == 2) {
                    filterUIs[[p]] <- selectInput(paste0("filter", i, "TSS"), "TSS", choices = c("off", "any", "either", "neither", "both"), selected = "off")
                } else {
                    filterUIs[[p]] <- selectInput(paste0("filter", i, "annote", j - 2), input$annotFiles$name[j - 2], choices = c("off", "any", "either", "neither", "both"), selected = "any")
                }
            }
        }
        tagList(filterUIs)
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

    maskTSSAnchors <- reactive({
        annotateAnchors(loops(), getHumanTSS(), "temp", maxgap = input$maxgap)@anchors |> mcols() |> _$temp
    })

    maskTSSLoopsLeft <- reactive({
        maskTSSAnchors()[loops()@interactions[, "left"]]
    })

    maskTSSLoopsRight <- reactive({
        maskTSSAnchors()[loops()@interactions[, "right"]]
    })

    maskTSSloopsList <- reactive({
        maskTSSloopsList <- vector('list', input$filterNum)
        for (i in seq_len(input$filterNum)) {
            mode <- input[[paste0("filter", i, "TSS")]]
            if (mode == "off") {
                maskTSSloopsList[[i]] <- rep(TRUE, nrow(loops()@rowData))
            } else if (mode == "any") {
                maskTSSloopsList[[i]] <-  maskTSSLoopsLeft() | maskTSSLoopsRight()
            } else if (mode == "either") {
                maskTSSloopsList[[i]] <- maskTSSLoopsLeft() & !maskTSSLoopsRight() | !maskTSSLoopsLeft() & maskTSSLoopsRight()
            } else if (mode == "neither") {
                maskTSSloopsList[[i]] <- !maskTSSLoopsLeft() & !maskTSSLoopsRight()
            } else if (mode == "both") {
                maskTSSloopsList[[i]] <- maskTSSLoopsLeft() & maskTSSLoopsRight()
            }
        }
        return(maskTSSloopsList)
    })

    annotGRangeList <- reactive({
        if (is.null(input$annotFiles)) {
            return(vector('list', 0))
        }
        annotGRangeList <- vector('list', nrow(input$annotFiles))
        for (i in seq_len(nrow(input$annotFiles))) {
            annotGRangeList[[i]] <- read_table(input$annotFiles$datapath[i], col_names = FALSE)[,1:3] |> setNames(c("chr", "start", "end")) |> makeGRangesFromDataFrame() |> rmchr()
        }
        return(annotGRangeList)
    })

    maskAnnotAnchorsList <- reactive({
        maskAnnotAnchorsList <- vector('list', length(annotGRangeList()))
        for (i in seq_len(length(annotGRangeList()))) {
            maskAnnotAnchorsList[[i]] <- annotateAnchors(loops(), annotGRangeList()[[i]], "temp", maxgap = input$maxgap)@anchors |> mcols() |> _$temp
        }
        return(maskAnnotAnchorsList)
    })

    maskAnnotLoopsLeft <- reactive({
        maskAnnotLoopsLeft <- vector('list', length(annotGRangeList()))
        for (i in seq_len(length(annotGRangeList()))) {
            maskAnnotLoopsLeft[[i]] <- maskAnnotAnchorsList()[[i]][loops()@interactions[, "left"]]
        }
        return(maskAnnotLoopsLeft)
    })

    maskAnnotLoopsRight <- reactive({
        maskAnnotLoopsRight <- vector('list', length(annotGRangeList()))
        for (i in seq_len(length(annotGRangeList()))) {
            maskAnnotLoopsRight[[i]] <- maskAnnotAnchorsList()[[i]][loops()@interactions[, "right"]]
        }
        return(maskAnnotLoopsRight)
    })

    maskAnnotLoopsFilterList <- reactive({
        maskAnnotLoopsFilterList <- vector('list', input$filterNum)
        for (i in seq_len(input$filterNum)) {
            maskAnnotLoopsFilterList[[i]] <- rep(TRUE, nrow(loops()@rowData))
            annoteNum <- ifelse(is.null(input$annotFiles), 0, nrow(input$annotFiles))
            for (j in seq_len(annoteNum)) {
                mode <- input[[paste0("filter", i, "annote", j)]]
                if (mode == "off") {
                    mask <- rep(TRUE, nrow(loops()@rowData))
                } else if (mode == "any") {
                    mask <- maskAnnotLoopsLeft()[[j]] | maskAnnotLoopsRight()[[j]]
                } else if (mode == "either") {
                    mask <- maskAnnotLoopsLeft()[[j]] & !maskAnnotLoopsRight()[[j]] | !maskAnnotLoopsLeft()[[j]] & maskAnnotLoopsRight()[[j]]
                } else if (mode == "neither") {
                    mask <- !maskAnnotLoopsLeft()[[j]] & !maskAnnotLoopsRight()[[j]]
                } else if (mode == "both") {
                    mask <- maskAnnotLoopsLeft()[[j]] & maskAnnotLoopsRight()[[j]]
                }
                maskAnnotLoopsFilterList[[i]] <- maskAnnotLoopsFilterList[[i]] & mask
            }
        }
        return(maskAnnotLoopsFilterList)
    })

    output$globalFilters <- renderUI({
        req(input$loopFiles)
        minLoopWidth <- min(loops()@rowData$loopWidth)
        maxLoopWidth <- max(loops()@rowData$loopWidth)
        widthFilter <- numericRangeInput("loopWidthRange", "loop width range", value = c(minLoopWidth, maxLoopWidth), min = minLoopWidth, max = maxLoopWidth, step = 1)
        maxLoopCount <- max(max(loops()@counts), 1)
        minLoopCount <- min(min(loops()@counts), maxLoopCount)
        countFilter <- numericRangeInput("loopCountRange", "loop count range", value = c(minLoopCount, maxLoopCount), min = minLoopCount, max = maxLoopCount, step = 1)
        tagList(widthFilter, countFilter)
    })

    maskLoopsWidth <- reactive({
        (loops()@rowData$loopWidth >= input$loopWidthRange[1]) & (loops()@rowData$loopWidth <= input$loopWidthRange[2])
    })

    maskLoopsCountList <- reactive({
        maskLoopsCountList <- vector('list', input$filterNum)
        for (i in seq_len(input$filterNum)) {
            if (length(input[[paste0("filter", i)]]) > 0) {
                selectedCounts <- loops()@counts[, input[[paste0("filter", i)]]]
                loopSelectedNum <- as.matrix(selectedCounts >= input$loopCountRange[1] & selectedCounts <= input$loopCountRange[2]) |> rowSums()
                maskLoopsCountList[[i]] <- loopSelectedNum > 0
            } else {
                maskLoopsCountList[[i]] <- rep(FALSE, nrow(loops()@rowData))
            }
        }
        return(maskLoopsCountList)
    })

    filterLoopList <- reactive({
        filterLoopList <- vector('list', input$filterNum)
        for (i in seq_len(input$filterNum)) {
            filterLoopList[[i]] <- loops() |> subsetLoops(maskFDRPValue() & !maskblackListLoops() & maskTSSloopsList()[[i]] & maskAnnotLoopsFilterList()[[i]] & maskLoopsWidth() & maskLoopsCountList()[[i]])
        }
        return(filterLoopList)
    })

    loopNums <- reactive({
        loopNums <- rep(0, input$filterNum)
        for (i in seq_len(input$filterNum)) {
            if (input$numberCount == "number") {
                loopNums[i] <- filterLoopList()[[i]]@counts[, input[[paste0("filter", i)]]] |> as.matrix() |> rowSums() |> as.logical() |> sum()
            } else if (input$numberCount == "count") {
                loopNums[i] <- filterLoopList()[[i]]@counts[, input[[paste0("filter", i)]]] |> sum()
            }
        }
        return(loopNums)
    })

    loopWidths <- reactive({
        loopWidthsList <- vector('list', input$filterNum)
        for (i in seq_len(input$filterNum)) {
            if (nrow(filterLoopList()[[i]]@rowData) == 0) {
                next
            }
            loopWidthsList[[i]] <- tibble(filter = paste0("filter", i), count = filterLoopList()[[i]]@counts[, input[[paste0("filter", i)]]] |> as.matrix() |> rowSums(), loopWidth = filterLoopList()[[i]]@rowData$loopWidth)
        }
        loopWidthsList |> bind_rows(tibble(filter = character(0), count = double(0), loopWidth = integer(0)))
    })

    hic <- reactive({
        HiCExperiment::import(input$hicFiles$datapath, focus = "chrM", resolution = resolutions() |> tail(n = 1))
    })

    reqInputs <- function() {
        req(input$loopFiles)
        for (i in seq_len(input$filterNum)) {
            req(input[[paste0("filter", i, "TSS")]])
            annoteNum <- ifelse(is.null(input$annotFiles), 0, nrow(input$annotFiles))
            for (j in seq_len(annoteNum)) {
                req(input[[paste0("filter", i, "annote", j)]])
            }
            req(input[[paste0("filter", i)]])
        }
    }

    loopNumGG <- reactive({
        tibble(filter = paste0("filter", seq_len(input$filterNum)), loopNum = loopNums()) |>
            ggplot(aes(filter, loopNum)) +
            geom_col()
    })

    output$loopNumPlot <- renderPlot({
        reqInputs()
        loopNumGG()
    })

    output$loopNumDownload <- downloadHandler(
        filename = "loopNum.pdf",
        content = function(file) {
            ggsave(file, loopNumGG())
        }
    )

    loopWidthGG <- reactive({
        loopWidths() |>
            ggplot(aes(x = loopWidth, weight = count, color = filter)) +
            stat_bin(geom = "line", bins = input$loopWidthBin, position = "identity")
    })

    output$loopWidthPlot <- renderPlot({
        reqInputs()
        loopWidthGG()
    })

    output$loopWidthDownload <- downloadHandler(
        filename = "loopWidth.pdf",
        content = function(file) {
            ggsave(file, loopWidthGG())
        }
    )

    loopWidthCdfGG <- reactive({
        loopWidths() |> 
            ggplot(aes(x = loopWidth, color = filter)) + 
            stat_ecdf(weight = count, pad = FALSE)
    })

    output$loopWidthCdfPlot <- renderPlot({
        reqInputs()
        loopWidthCdfGG()
    })

    output$loopWidthCdfDownload <- downloadHandler(
        filename = "loopWidthCdf.pdf",
        content = function(file) {
            ggsave(file, loopWidthCdfGG())
        }
    )

    bedpeFiles <- reactive({
        bedpeFiles <- vector('list', input$filterNum)
        for (i in seq_len(input$filterNum)) {
            loops <- filterLoopList()[[i]]
            seqnames <- loops@anchors |> rmchr() |> addchr() |> _@seqnames
            chrs <- rep(levels(seqnames@values), times = c(seqnames@lengths, rep(0, length(levels(seqnames@values)) - length(seqnames@lengths))))
            anchors <- loops@anchors |> rmchr() |> addchr() |> as_tibble() |> _[, 1:3]
            leftAnchors <- anchors |> _[loops@interactions[, 'left'], ] |> setNames(c("chrom1", "start1", "end1"))
            rightAnchors <- anchors |> _[loops@interactions[, 'right'], ] |> setNames(c("chrom2", "start2", "end2"))
            bedpeFiles[[i]] <- tempfile(fileext = ".bedpe")
            cbind(leftAnchors, rightAnchors) |> mutate(name = ".", count = loops@counts[, input[[paste0("filter", i)]]] |> as.matrix() |> rowSums()) |> write_tsv(file = bedpeFiles[[i]], col_names = FALSE)
        }
        return(bedpeFiles)
    })

    output$heatMapPlot <- renderUI({
        req(input$hicFiles)
        plotOutputList <- vector('list', 2 * input$filterNum)
        for (i in seq_len(input$filterNum)) {
            plotOutputList[[2 * i - 1]] <- plotOutput(paste0("heatMap", i), width = "1000px", height = "1000px")
            plotOutputList[[2 * i]] <- downloadButton(paste0("heatMap", i, "download"))
        }
        tagList(plotOutputList)
    })

    heatMapGGList <- reactive({
        balanced <- ifelse(input$balanced, "balanced", "count")
        maxDistance <- switch(input$horizontal + 1, NULL, input$maxDistance)
        heatMapGGList <- vector('list', input$filterNum)
        for (i in seq_len(input$filterNum)) {
            loops <- bedpeFiles()[[i]] |> import() |> InteractionSet::makeGInteractionsFromGRangesPairs()
            heatMapGGList[[i]] <- hic() |> refocus(paste0(input$chromosome, ":", input$start, "-", input$end)) |> zoom(as.integer(input$resolution)) |> plotMatrix(use.scores = balanced, scale = input$scale, loops = loops, maxDistance = maxDistance)
        }
        return(heatMapGGList)
    })

    observe({
        reqInputs()
        req(input$hicFiles)
        req(input$chromosome)
        req(input$start)
        req(input$end)
        req(input$resolution)
        for (i in seq_len(input$filterNum)) {
            local({
                my_i <- i
                output[[paste0("heatMap", my_i)]] <- renderPlot({
                    isolate({
                        heatMapGGList()[[my_i]]
                    })
                })
                output[[paste0("heatMap", my_i, "download")]] <- downloadHandler(
                    filename = paste0("heatMap", my_i, ".pdf"),
                    content = function(file) {
                        ggsave(file, heatMapGGList()[[my_i]])
                    }
                )
            })
        }
    }) |> bindEvent(input$renderHeat)

    output$aggregatePlot <- renderUI({
        req(input$hicFiles)
        imageOutputList = vector("list", 2 * input$filterNum)
        for (i in seq_len(input$filterNum)) {
            imageOutputList[[2 * i - 1]] <- imageOutput(paste0("aggregate", i), width = "400px", height = "400px")
            imageOutputList[[2 * i]] <- downloadButton(paste0("aggregate", i, "download"))
        }
        tagList(imageOutputList)
    })

    observe({
        reqInputs()
        req(input$hicFiles)
        req(input$aggResolution)
        flankBps <- input$flank * as.integer(input$aggResolution)
        for (i in seq_len(input$filterNum)) {
            local({
                my_i <- i
                clpyFile <- tempfile(fileext = ".clpy")
                pngFile <- sub(".clpy$", ".png", clpyFile)
                if (!is.null(input$view)) {
                    coolpupCmd = sprintf("coolpup.py --view %s --flank %d -o %s --seed 0 --weight-name %s %s::/resolutions/%s %s", input$view$datapath, flankBps, clpyFile, input$normalization, input$hicFiles$datapath, input$aggResolution, bedpeFiles()[[my_i]])
                } else {
                    coolpupCmd = sprintf("coolpup.py --flank %d -o %s --seed 0 --weight-name %s %s::/resolutions/%s %s", flankBps, clpyFile, input$normalization, input$hicFiles$datapath, input$aggResolution, bedpeFiles()[[my_i]])
                }
                coolpupCmd |> system()
                if (input$setVminVmax) {
                    req(input$vmin)
                    req(input$vmax)
                    plotpupCmd = sprintf("plotpup.py --input_pups %s -o %s --vmin %f --vmax %f", clpyFile, pngFile, input$vmin, input$vmax)
                } else {
                    plotpupCmd = sprintf("plotpup.py --input_pups %s -o %s", clpyFile, pngFile)
                }
                plotpupCmd |> system()
                output[[paste0("aggregate", my_i)]] <- renderImage(
                    {
                        list(src = pngFile)
                    },
                    deleteFile = FALSE 
                )
                output[[paste0("aggregate", my_i, "download")]] <- downloadHandler(
                    filename = paste0("aggregate", my_i, ".pdf"),
                    content = function(file) {
                        if (input$setVminVmax) {
                            req(input$vmin)
                            req(input$vmax)
                            plotpupCmd = sprintf("plotpup.py --input_pups %s -o %s --vmin %f --vmax %f", clpyFile, file, input$vmin, input$vmax)
                        } else {
                            plotpupCmd = sprintf("plotpup.py --input_pups %s -o %s", clpyFile, file)
                        }
                        plotpupCmd |> system()
                    }
                )
            })
        }
    }) |> bindEvent(input$renderAggregate)
}

shinyApp(ui = ui, server = server)
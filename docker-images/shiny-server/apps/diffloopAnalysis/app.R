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
library(halfmoon)

chromInfo <- getChromInfoFromUCSC("hg19")

source("../helpers/mangoFDRPValue.R")

library(spatstat)

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
        uiOutput("loopNumPlot")
    ),
    tabPanel("loop width",
        numericInput("loopWidthBin", "loop width bin number", value = 20, min = 0, step = 1),
        uiOutput("loopWidthPlot")
    ),
    tabPanel("loop width cumulative distribution",
        uiOutput("loopWidthCdfPlot")
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

    resolutions <- reactive({
        resolutions <- availableResolutions(input$hicFiles$datapath)
        return(resolutions[resolutions >= 10000])
    })

    observe({
        toggleState("vmin")
        toggleState("vmax")
    }) |> bindEvent(input$setVminVmax)

    # Use proxy$xxxxx as a proxy of input$xxxxx. update??????Input will not update input$xxxxx until the client send the updated values back to the server. proxy$xxxxx helps to block renderUI until input$xxxxx got updated.
    proxy <- reactiveValues(
        resolution = NULL,
        aggResolution = NULL,
        chromosome = NULL,
        start = NULL,
        end = NULL
    )
    observe({
        for (name in names(proxy)) {
            proxy[[name]] <- input[[name]]
        }
    })
    observe({
        req(input$hicFiles)
        proxy$resolution <- NULL
        proxy$aggResolution <- NULL
        proxy$chromosome <- NULL
        updateSelectInput(inputId = "resolution", choices = resolutions())
        updateSelectInput(inputId = "aggResolution", choices = resolutions())
        updateSelectInput(inputId = "chromosome", choices = availableChromosomes(input$hicFiles$datapath)@seqnames)
    })

    observe({
        req(proxy$chromosome)
        proxy$start <- NULL
        proxy$end <- NULL
        chromSize <- chromInfo$size[chromInfo$chrom == proxy$chromosome]
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
        maxLoopCount <- max(max(loops()@counts), 2)
        minLoopCount <- min(min(loops()@counts), maxLoopCount)
        countFilter <- numericRangeInput("loopCountRange", "loop count range", value = c(max(minLoopCount, 2), maxLoopCount), min = minLoopCount, max = maxLoopCount, step = 1)
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

    loopNumTempFile <- tempfile(tmpdir=file.path("www", session$token))
    output$loopNumPlot <- renderUI({
        reqInputs()
        ggObj <- tibble(filter = paste0("filter", seq_len(input$filterNum)), loopNum = loopNums()) |>
            ggplot(aes(filter, loopNum)) +
            geom_col()
        ggsave(paste0(loopNumTempFile, ".pdf"), plot = ggObj)
        tags$iframe(src = paste0(sub("^www/", "", loopNumTempFile), ".pdf"), height = "1200px", width = "100%")
    })

    loopWidthTempFile <- tempfile(tmpdir=file.path("www", session$token))
    output$loopWidthPlot <- renderUI({
        reqInputs()
        ggObj <- loopWidths() |>
            ggplot(aes(x = loopWidth, weight = count, color = filter)) +
            stat_bin(geom = "line", bins = input$loopWidthBin, position = "identity")
        ggsave(paste0(loopWidthTempFile, ".pdf"), plot = ggObj)
        tags$iframe(src = paste0(sub("^www/", "", loopWidthTempFile), ".pdf"), height = "1200px", width = "100%")
    })

    loopWidthCdfTempFile <- tempfile(tmpdir=file.path("www", session$token))
    output$loopWidthCdfPlot <- renderUI({
        reqInputs()
        ggObj <- loopWidths() |> 
            ggplot(aes(x = loopWidth, color = filter)) +
            geom_ecdf(aes(weights = count))
        ggsave(paste0(loopWidthCdfTempFile, ".pdf"), plot = ggObj)
        tags$iframe(src = paste0(sub("^www/", "", loopWidthCdfTempFile), ".pdf"), height = "1200px", width = "100%")
    })

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

    heatMapGGList <- reactive({
        balanced <- ifelse(input$balanced, "balanced", "count")
        maxDistance <- switch(input$horizontal + 1, NULL, input$maxDistance)
        heatMapGGList <- vector('list', input$filterNum)
        for (i in seq_len(input$filterNum)) {
            loops <- bedpeFiles()[[i]] |> import() |> InteractionSet::makeGInteractionsFromGRangesPairs()
            heatMapGGList[[i]] <- hic() |> refocus(paste0(proxy$chromosome, ":", proxy$start, "-", proxy$end)) |> zoom(as.integer(proxy$resolution)) |> plotMatrix(use.scores = balanced, scale = input$scale, loops = loops, maxDistance = maxDistance)
        }
        return(heatMapGGList)
    })

    heatMapTempFile <- tempfile(tmpdir=file.path("www", session$token))
    observe({
        reqInputs()
        req(input$hicFiles)
        req(proxy$chromosome)
        req(proxy$start)
        req(proxy$end)
        req(proxy$resolution)

        # render UIs for holding iframe
        output$heatMapPlot <- renderUI({
            # isolate prevents renderUI from automatic update, it is update only if bindEvent is satisfied
            isolate({
                plotOutputList <- vector('list', input$filterNum)
                for (i in seq_len(input$filterNum)) {
                    plotOutputList[[i]] <- uiOutput(paste0("heatMap", i))
                }
                tagList(plotOutputList)
            })
        })

        # core code to generate iframe
        for (i in seq_len(input$filterNum)) {
            local({
                my_i <- i
                output[[paste0("heatMap", my_i)]] <- renderUI({
                    # isolate prevents renderUI from automatic update, it is update only if bindEvent is satisfied
                    isolate({
                        heatMapTempFileIter <- paste0(heatMapTempFile, ".", my_i, ".pdf")
                        ggsave(heatMapTempFileIter, plot = heatMapGGList()[[my_i]])
                        tags$iframe(src = sub("^www/", "", heatMapTempFileIter), height = "1200px", width = "100%")
                    })
                })
            })
        }
    }) |> bindEvent(input$renderHeat)

    aggregateTempFile <- tempfile(tmpdir=file.path("www", session$token))
    observe({
        reqInputs()
        req(input$hicFiles)
        req(proxy$aggResolution)

        # render UIs for holding iframe
        output$aggregatePlot <- renderUI({
            # isolate prevents renderUI from automatic update, it is update only if bindEvent is satisfied
            isolate({
                imageOutputList = vector("list", input$filterNum)
                for (i in seq_len(input$filterNum)) {
                    imageOutputList[[i]] <- uiOutput(paste0("aggregate", i))
                }
                tagList(imageOutputList)
            })
        })

        # core code to generate iframe
        flankBps <- input$flank * as.integer(proxy$aggResolution)
        for (i in seq_len(input$filterNum)) {
            local({
                my_i <- i
                clpyFile <- paste0(aggregateTempFile, ".", my_i, ".clpy")
                pdfFile <- sub(".clpy$", ".pdf", clpyFile)
                if (!is.null(input$view)) {
                    coolpupCmd = sprintf("coolpup.py --view %s --flank %d -o %s --seed 0 --weight-name %s %s::/resolutions/%s %s", input$view$datapath, flankBps, clpyFile, input$normalization, input$hicFiles$datapath, proxy$aggResolution, bedpeFiles()[[my_i]])
                } else {
                    coolpupCmd = sprintf("coolpup.py --flank %d -o %s --seed 0 --weight-name %s %s::/resolutions/%s %s", flankBps, clpyFile, input$normalization, input$hicFiles$datapath, proxy$aggResolution, bedpeFiles()[[my_i]])
                }
                coolpupCmd |> system()
                if (input$setVminVmax) {
                    req(input$vmin)
                    req(input$vmax)
                    plotpupCmd = sprintf("plotpup.py --input_pups %s -o %s --vmin %f --vmax %f", clpyFile, pdfFile, input$vmin, input$vmax)
                } else {
                    plotpupCmd = sprintf("plotpup.py --input_pups %s -o %s", clpyFile, pdfFile)
                }
                plotpupCmd |> system()
                output[[paste0("aggregate", my_i)]] <- renderUI({
                    # isolate prevents renderUI from automatic update, it is update only if bindEvent is satisfied
                    isolate({
                        tags$iframe(src = sub("^www/", "", pdfFile), height = "1200px", width = "100%")
                    })
                })
            })
        }
    }) |> bindEvent(input$renderAggregate)
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

shinyApp(ui = ui, server = server)
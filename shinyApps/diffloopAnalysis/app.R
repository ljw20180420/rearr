library(shiny)
library(bslib)
library(tidyverse)
library(diffloop)
library(GenomicRanges)

source("../helpers/mangoFDRPValue.R")
source("../helpers/stat_ecdf_weighted.R")

options(shiny.maxRequestSize = 1000 * 1024^2)

ui <- page_sidebar(
    title = "Diffloop Analysis",
    sidebar = sidebar(
        helpText("This App use diffloop to analyze loop data."),
        fileInput("loopFiles", label = "Loop files (bedpe or mango)", multiple = TRUE),
        numericInput("mergegap", label = "merge gap", value = 500),
        numericInput("nbins", label = "nbins", value = 10, min = 0, step = 1),

        numericInput("FDR", label = "FDR", value = 1, min = 0, max = 1, step = 0.01),
        numericInput("PValue", label = "PValue", value = 1, min = 0, max = 1, step = 0.01),

        numericInput("maxgap", label = "max gap between blackList/annotation and anchor", value = 0, min = 0, step = 100),

        fileInput("blackListFiles", label = "black list regions (bed files)", multiple = TRUE),

        fileInput("annotFiles", label = "annotation files", multiple = TRUE),

        sliderInput("loopWidthRange", label = "loop width range", min = 0, max = 0, value = c(0, 0), step = 1),

        sliderInput("loopCountRange", label = "loop count range", min = 1, max = 1, value = c(1, 1), step = 1),

        fileInput("hicFiles", label = "hic files", multiple = TRUE),

        actionButton("addFilter", "add filter"),
        actionButton("clearFilter", "clear filter"),

        selectInput("outputType", "output type", choices = c("loopNum", "loopWidth", "loopWidthCdf", "heatMap"), selected = "loopNum")
    ),
    plotOutput("anyPlot")
)

server <- function(input, output) {
    beddir <- tempdir()
    filterNum <- reactiveVal(0)
    annoteNum <- reactiveVal(0)

    samples <- reactive({
        input$loopFiles$name |> strsplit('.', fixed = TRUE) |> vapply(function(x) x[1], "")
    })

    observeEvent(input$addFilter, {
        insertUI(
            selector = "#addFilter",
            where = "afterEnd",
            ui = selectInput(
                inputId = paste0("filter", filterNum() + 1, "TSS"),
                label = "TSS",
                choices = c("off", "any", "either", "neither", "both"),
                selected = "any"
            ),
            immediate = TRUE
        )
        for (i in seq_len(annoteNum())) {
            insertUI(
                selector = "#addFilter",
                where = "afterEnd",
                ui = selectInput(
                    inputId = paste0("filter", filterNum() + 1, "annote", i),
                    label = input$annotFiles$name[i],
                    choices = c("off", "any", "either", "neither", "both"),
                    selected = "any"
                ),
                immediate = TRUE
            )
        }
        insertUI(
            selector = "#addFilter",
            where = "afterEnd",
            ui = checkboxGroupInput(
                inputId = paste0("filter", filterNum() + 1),
                label = paste0("filter", filterNum() + 1),
                choices = samples(),
                selected = samples()
            ),
            immediate = TRUE
        )
        filterNum(filterNum() + 1)
    })

    observeEvent(ignoreInit = TRUE, list(
            input$loopFiles,
            input$annotFiles,
            input$clearFilter
        ),
        {
            for (i in seq_len(filterNum())) {
                removeUI(selector = sprintf("div:has(>> #%s)", paste0("filter", i, "TSS")), immediate = TRUE)
                for (j in seq_len(annoteNum())) {
                    removeUI(selector = sprintf("div:has(>> #%s)", paste0("filter", i, "annote", j)), immediate = TRUE)
                }
                removeUI(selector = paste0("#filter", i), immediate = TRUE)
            }
            if (is.null(input$annotFiles)) {
                annoteNum(0)
            } else {
                annoteNum(nrow(input$annotFiles))
            }
            filterNum(0)
        }
    )

    observeEvent(input$outputType, {
        if (input$outputType == "loopWidth") {
            insertUI(
                selector = "div:has(>> #outputType)",
                where = "afterEnd",
                ui = numericInput(
                    inputId = "loopWidthBin",
                    label = "loop width bin number",
                    value = 20,
                    min = 0,
                    step = 1
                ),
                immediate = TRUE
            )
        } 
        
        if (input$outputType != "loopWidth") {
            removeUI(selector = sprintf('.shiny-input-container:has(#%s)','loopWidthBin'), immediate = TRUE)
        }

        if (input$outputType == "heatMap" && !is.null(input$hicFiles)) {
            insertUI(
                selector = "div:has(>> #outputType)",
                where = "afterEnd",
                ui = selectInput(
                    inputId = "resolution",
                    label = "resolution",
                    choices = availableResolutions(input$hicFiles$datapath[1])
                ),
                immediate = TRUE
            )
            insertUI(
                selector = "div:has(>> #outputType)",
                where = "afterEnd",
                ui = checkboxInput(
                    inputId = "balanced",
                    label = "balanced",
                    value = TRUE
                ),
                immediate = TRUE
            )
            insertUI(
                selector = "div:has(>> #outputType)",
                where = "afterEnd",
                ui = numericInput(
                    inputId = "end",
                    label = "end",
                    value = NA,
                    min = 0
                ),
                immediate = TRUE
            )
            insertUI(
                selector = "div:has(>> #outputType)",
                where = "afterEnd",
                ui = numericInput(
                    inputId = "start",
                    label = "start",
                    value = NA,
                    min = 0
                ),
                immediate = TRUE
            )
            insertUI(
                selector = "div:has(>> #outputType)",
                where = "afterEnd",
                ui = textInput(
                    inputId = "chromosome",
                    label = "chromosome"
                ),
                immediate = TRUE
            )
        }

        if (input$outputType != "heatMap") {
            removeUI(selector = sprintf('.shiny-input-container:has(#%s)', "resolution"), immediate = TRUE)
            removeUI(selector = sprintf('.shiny-input-container:has(#%s)', "balanced"), immediate = TRUE)
            removeUI(selector = sprintf('.shiny-input-container:has(#%s)', "end"), immediate = TRUE)
            removeUI(selector = sprintf('.shiny-input-container:has(#%s)', "start"), immediate = TRUE)
            removeUI(selector = sprintf('.shiny-input-container:has(#%s)', "chromosome"), immediate = TRUE)
        }
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
        maskTSSloopsList <- vector('list', filterNum())
        for (i in seq_len(filterNum())) {
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
        maskAnnotLoopsFilterList <- vector('list', filterNum())
        for (i in seq_len(filterNum())) {
            maskAnnotLoopsFilterList[[i]] <- rep(TRUE, nrow(loops()@rowData))
            for (j in seq_len(annoteNum())) {
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

    maskLoopsCountList <- reactive({
        maskLoopsCountList <- vector('list', filterNum())
        for (i in seq_len(filterNum())) {
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
        filterLoopList <- vector('list', filterNum())
        for (i in seq_len(filterNum())) {
            filterLoopList[[i]] <- loops() |> subsetLoops(maskFDRPValue() & !maskblackListLoops() & maskTSSloopsList()[[i]] & maskAnnotLoopsFilterList()[[i]] & maskLoopsWidth() & maskLoopsCountList()[[i]])
        }
        return(filterLoopList)
    })

    loopNums <- reactive({
        loopNums <- rep(0, filterNum())
        for (i in seq_len(filterNum())) {
            loopNums[i] <- filterLoopList()[[i]]@counts[, input[[paste0("filter", i)]]] |> sum()
        }
        return(loopNums)
    })

    loopWidths <- reactive({
        loopWidthsList <- vector('list', filterNum())
        for (i in seq_len(filterNum())) {
            if (nrow(filterLoopList()[[i]]@rowData) == 0) {
                next
            }
            loopWidthsList[[i]] <- tibble(filter = paste0("filter", i), count = filterLoopList()[[i]]@counts |> rowSums(), loopWidth = filterLoopList()[[i]]@rowData$loopWidth)
        }
        loopWidthsList |> bind_rows(tibble(filter = character(0), count = double(0), loopWidth = integer(0)))
    })

    hics <- reactive({
        hics <- vector('list', nrow(input$hicFiles))
        for (i in seq_len(nrow(input$hicFiles))) {
            hics[[i]] <- HiCExperiment::import(input$hicFiles$datapath[i],
                focus = paste0(input$chromosome, ":", input$start, "-", input$end),
                resolution = as.integer(input$resolution)
            )
        }
    })

    output$anyPlot <- renderPlot({
        req(filterNum() > 0)
        req(input$loopFiles)
        for (i in seq_len(filterNum())) {
            req(input[[paste0("filter", i, "TSS")]])
            for (j in seq_len(annoteNum())) {
                req(input[[paste0("filter", i, "annote", j)]])
            }
            req(input[[paste0("filter", i)]])
        }

        if (input$outputType == "loopNum") {
            tibble(filter = paste0("filter", seq_len(filterNum())), loopNum = loopNums()) |>
                ggplot(aes(filter, loopNum)) +
                geom_col()
        } else if (input$outputType == "loopWidth") {
            req(input$loopWidthBin)
            loopWidths() |>
                ggplot(aes(x = loopWidth, weight = count, color = filter)) +
                stat_bin(geom = "line", bins = input$loopWidthBin, position = "identity")
        } else if (input$outputType == "loopWidthCdf") {
            loopWidths() |> 
                ggplot(aes(x = loopWidth, color = filter)) + 
                stat_ecdf(weight = count, pad = FALSE)
        } else if (input$outputType == "heatMap") {
            req(input$hicFiles)
            req(input$chromosome)
            req(input$start)
            req(input$end)
            req(input$resolution)
            req(input$balanced)
            balanced <- ifelse(input$balanced, "balanced", "count")
            plotMatrix(hics()[[1]])
        }
    })
}

shinyApp(ui = ui, server = server)
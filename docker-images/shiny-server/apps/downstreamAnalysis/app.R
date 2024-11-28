library(shiny)
library(bslib)
library(tidyverse)
library(shinyWidgets)
library(reshape2)
library(scales)
library(ggseqlogo)
library(ggforce)
library(waffle)

options(shiny.maxRequestSize = 100 * 1024^3)

dir(path = "helpers", full.names = TRUE) |> lapply(source)

# Define UI ----
ui <- navbarPage(
    "downstream analysis",
    theme = bslib::bs_theme(version = 5),
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
        fileInput("sgRNAfile", "sgRNA file"),
        selectInput("editTarget", "edit target", choices = c("templated insertion", "random insertion", "insertion", "deletion", "templated indel", "indel", "wild type"))
    ),
    tabPanel("browser",
        uiOutput("browserReadRangeUI"),
        htmlOutput("alignments", class = "alignments", inline = TRUE)
    ),
    tabPanel("one",
        textInput("alignOneRef1", label = "reference1", value = ""),
        textInput("alignOneRef2", label = "reference2", value = ""),
        numericInput("alignOneCut1", "cut1", value = NA, min = 0),
        numericInput("alignOneCut2", "cut2", value = NA, min = 0),
        selectInput("alignOnePAM1", "PAM1", choices = c("NGG", "CCN"), selected = "NGG"),
        selectInput("alignOnePAM2", "PAM2", choices = c("NGG", "CCN"), selected = "NGG"),
        textInput("alignOneQuery", label = "query read", value = ""),
        htmlOutput("alignPair", class = "alignments", inline = TRUE)
    ),
    tabPanel("base",
        uiOutput("baseSubFreqPlot")
    ),
    tabPanel("position",
        selectInput("positionalMode", "display mode", choices = c("histgram base", "histgram indel", "read base", "read snp", "logo probability", "logo bits", "logo custom"), selected = "histgram base"),
        uiOutput("posBaseRef1Plot"),
        uiOutput("posBaseRef2Plot")
    ),
    tabPanel("MMEJ",
        selectInput("microRefId", "reference id", choices = NULL),
        selectInput("microMode", "microhomology count display mode", choices = c("repeat", "separate")),
        numericInput("microThres", "minimal micro homology length", value = 4, min = 1),
        uiOutput("mhMatrixPlot")
    ),
    tabPanel("classify",
        checkboxInput("claClaDistinctTemp", "discriminate templated insertion and random insertion"),
        selectInput("claClaMode", "mode", choices = c("pie", "waffle")),
        uiOutput("claClaPlot")
    ),
    tabPanel("distribution",
        selectInput("distriTarget", "what", choices = c("score", "cut1", "cut2", "ref1Len", "ref2Len", "ref1End", "ref2Start", "upInsert", "downInsert", "randInsert", "templatedInsert", "insert", "upDelete", "downDelete", "delete"), multiple = TRUE),
        selectInput("distriMode", "mode", choices = c("continuous", "discrete")),
        uiOutput("distriPlot")
    ),
    tabPanel("pairwise",
        selectInput("pairwiseX", "x", choices = c("score", "cut1", "cut2", "ref1Len", "ref2Len", "ref1End", "ref2Start", "upInsert", "downInsert", "randInsert", "templatedInsert", "insert", "upDelete", "downDelete", "delete")),
        selectInput("pairwiseY", "y", choices = c("score", "cut1", "cut2", "ref1Len", "ref2Len", "ref1End", "ref2Start", "upInsert", "downInsert", "randInsert", "templatedInsert", "insert", "upDelete", "downDelete", "delete")),
        selectInput("pairwiseXscale", "x scale", choices = c("identity", "log10")),
        selectInput("pairwiseYscale", "y scale", choices = c("identity", "log10")),
        selectInput("pairwiseMethod", "method", choices = c("auto", "lm", "glm", "gam", "loess")),
        numericInput("pairwiseSpan", "span", 0.75, min = 0, max = 1, step = 0.01),
        uiOutput("pairwisePlot"),
        textOutput("pairwiseWarning")
    ),
    tabPanel("polygon",
        uiOutput("polyInsert1Plot"),
        uiOutput("polyInsert2Plot")
    ),
    tabPanel("arc",
        uiOutput("arcDelete1Plot"),
        uiOutput("arcDelete2Plot")
    ),
    tabPanel("kpLogo",
        numericInput("kpLogoKmer", "kmer", 1, min = 1, max = 10, step = 1),
        sliderInput("kpLogoRegion", "region", 0, 0, c(0, 0)),
        selectInput("kpLogoMethod", "method", choices = c("weight", "background")),
        numericInput("kpLogoCountThres", "count threshold", 0, min = 0),
        uiOutput("kpLogoIframe")
    ),
    tabPanel("kmer",
        uiOutput("kmerRangeUI"),
        uiOutput("kmerIframe")
    )
)

# Define server logic ----
server <- function(input, output, session) {
    ################################
    # session start
    ################################
    file.path("www", session$token) |> dir.create()
    ################################
    # session end
    ################################
    onSessionEnded(function() {
        unlink(file.path("www", session$token), recursive = TRUE)
    })
    ################################
    # sidebar
    ################################
    sgRNAs <- reactive({
        readLines(input$sgRNAfile$datapath)
    })
    algTibble <- reactive({
        algLines <- lapply(
            input$algfiles$datapath,
            function(datapath) {
                readLines(gzfile(datapath))
                }
            ) |> unlist()
        algLines[seq(1, length(algLines), 3)] |> I() |> read_tsv(
            col_names = c("index", "count", "score", "refId", "upDangle", "ref1Start", "query1Start", "ref1End", "query1End", "randInsert", "ref2Start", "query2Start", "ref2End", "query2End", "downDangle", "cut1", "cut2"),
            col_types = "iidiciiiiciiiicii",
            col_select = -1,
            na = character()
        ) |>
        mutate(
            refLine = algLines[seq(2, length(algLines), 3)],
            queryLine = algLines[seq(3, length(algLines), 3)],
            refNoGap = refLine |> str_replace_all("-", ""),
            ref1Len = refNoGap |> vapply(function(x) gregexpr("[acgtn]", x)[[1]][2], 0, USE.NAMES = FALSE),
            cut2 = cut2 - ref1Len,
            ref2Start = ref2Start - ref1Len,
            ref2End = ref2End - ref1Len,
            ref2Len = nchar(refNoGap) - ref1Len,
            upInsert = ifelse(ref1End > cut1, ref1End - cut1, 0),
            upDelete = ifelse(ref1End < cut1, cut1 - ref1End, 0),
            downInsert = ifelse(ref2Start < cut2, cut2 - ref2Start, 0),
            downDelete = ifelse(ref2Start > cut2, ref2Start - cut2, 0),
            templatedInsert = upInsert + downInsert,
            delete = upDelete + downDelete,
            insert = nchar(randInsert) + templatedInsert
        )
    })
    algMetaData <- reactive({
        tibble(
            totalCount = sum(algTibble()$count),
            maxCut1 = max(algTibble()$cut1),
            maxCut2 = max(algTibble()$cut2),
            maxCut1down = max(algTibble()$ref1Len - algTibble()$cut1),
            maxCut2down = max(algTibble()$ref2Len - algTibble()$cut2),
            maxRandInsert = max(nchar(algTibble()$randInsert))
        )
    })
    editTarget <- reactive({
        if (input$editTarget == "templated insertion") {
            return(as.logical(algTibble()$templatedInsert))
        } else if (input$editTarget == "random insertion") {
            return(as.logical(algTibble()$randInsert))
        } else if (input$editTarget == "insertion") {
            return(as.logical(algTibble()$insert))
        } else if (input$editTarget == "deletion") {
            return(as.logical(algTibble()$delete))
        } else if (input$editTarget == "templated indel") {
            return(algTibble()$templatedInsert & algTibble()$delete)
        } else if (input$editTarget == "indel") {
            return(algTibble()$insert & algTibble()$delete)
        } else if (input$editTarget == "wild type") {
            return(!algTibble()$insert & !algTibble()$delete)
        }
    })

    ################################
    # alignment browser
    ################################
    # use browserReadRange as a proxy of input$browserReadRange to prevent render output$alignments twice (one before input$browserReadRange is returned by client, and one after that)
    browserReadRange <- reactiveVal(NULL)
    observe({
        browserReadRange(input$browserReadRange)
    })
    output$browserReadRangeUI <- renderUI({
        req(input$algfiles)
        browserReadRange(NULL)
        numericRangeInput("browserReadRange", "read range to browse", value = c(1, min(1000, nrow(algTibble()))), min = 1, max = nrow(algTibble()), step = 1)
    })
    output$alignments <- renderText({
        req(input$algfiles)
        req(browserReadRange())
        paste(getMarkdownFromAlign(algTibble()[browserReadRange()[1]:browserReadRange()[2],]), collapse = "")
    })

    #########################################################
    # single alignment
    #########################################################
    algOneRefFile <- tempfile(tmpdir=file.path("www", session$token))
    algOneAlgFile <- tempfile(tmpdir=file.path("www", session$token))
    output$alignPair <- renderText({
        req(input$alignOneRef1)
        req(input$alignOneRef2)
        req(input$alignOneCut1)
        req(input$alignOneCut2)
        req(input$alignOneQuery)
        alignOne(input$alignOneRef1, input$alignOneRef2, input$alignOneCut1, input$alignOneCut2, input$alignOnePAM1, input$alignOnePAM2, input$alignOneQuery, algOneRefFile, algOneAlgFile)
    })

    #########################################################
    # base substitution frequency
    #########################################################
    baseSubFreqTempFile <- tempfile(tmpdir=file.path("www", session$token), fileext=".pdf")
    output$baseSubFreqPlot <- renderUI({
        req(input$algfiles)
        ggFig <- algTibble() |> countBaseSubstitute() |>
            ggplot(aes(sub, count)) +
            geom_col() +
            scale_y_continuous(expand = c(0, 0))
        ggsave(baseSubFreqTempFile, plot = ggFig)
        tags$iframe(src = sub("^www/", "", baseSubFreqTempFile), height = "1200px", width = "100%")
    })

    #################################
    # positional statistics
    #################################
    refGapMat <- reactive({
        algTibble()$refLine |> extendToSameLength()
    })
    ref1Mat <- reactive({
        algTibble()$refNoGap |> substr(1, algTibble()$ref1Len) |> extendToAlignCut(algTibble()$cut1)
    })
    ref2Mat <- reactive({
        algTibble()$refNoGap |> substr(algTibble()$ref1Len + 1, algTibble()$ref1Len + algTibble()$ref2Len) |> extendToAlignCut(algTibble()$cut2)
    })
    queryAllSeqs <- reactive({
        tmpMat <- algTibble()$queryLine |> extendToSameLength() |> t() |> as.vector()
        vectorToStringVector(tmpMat[(refGapMat() |> t() |> as.vector()) != '-'], algTibble()$ref1Len + algTibble()$ref2Len)
    })
    query1Mat <- reactive({
        queryAllSeqs() |> substr(1, algTibble()$ref1Len) |> extendToAlignCut(algTibble()$cut1)
    })
    query2Mat <- reactive({
        queryAllSeqs() |> substr(algTibble()$ref1Len + 1, nchar(queryAllSeqs())) |> extendToAlignCut(algTibble()$cut2)
    })
    base1Freq <- reactive({
        getPositionalBaseFreq(query1Mat(), algTibble()$count)
    })
    base2Freq <- reactive({
        getPositionalBaseFreq(query2Mat(), algTibble()$count)
    })
    base1Tibble <- reactive({
        base1Freq() |> posMatrixToTibble(algMetaData()$maxCut1)
    })
    base2Tibble <- reactive({
        base2Freq() |> posMatrixToTibble(algMetaData()$maxCut2)
    })
    MSD1Tibble <- reactive({
        getPositionalMSDTibble(query1Mat(), ref1Mat(), algTibble()$count, algMetaData()$maxCut1)
    })
    MSD2Tibble <- reactive({
        getPositionalMSDTibble(query2Mat(), ref2Mat(), algTibble()$count, algMetaData()$maxCut2)
    })
    insert1Count <- reactive({
        algTibble()$refLine |> vapply(function(x) substr(x, 1, gregexpr("[acgtn]", x)[[1]][3] - 1), "", USE.NAMES = FALSE) |> strsplit("") |> calInsertionCount(cuts = algTibble()$cut1, maxCutDown = algMetaData()$maxCut1down)
    })
    insert2Count <- reactive({
        algTibble()$refLine |> vapply(function(x) substr(x, gregexpr("[acgtn]", x)[[1]][2] + 1, nchar(x)), "", USE.NAMES = FALSE) |> strsplit("") |> calInsertionCount(cuts = algTibble()$cut2, maxCutDown = algMetaData()$maxCut2down)
    })
    read1Tibble <- reactive({
        getPositionalReads(query1Mat(), algTibble()$count, algMetaData()$totalCount * 0.001, algMetaData()$maxCut1)
    })
    read2Tibble <- reactive({
        getPositionalReads(query2Mat(), algTibble()$count, algMetaData()$totalCount * 0.001, algMetaData()$maxCut2)
    })
    snp1Tibble <- reactive({
        queryMat <- query1Mat()
        refMat <- toupper(ref1Mat())
        queryMat[queryMat != refMat & queryMat != "-"] = "S"
        queryMat[queryMat == refMat] = "M"
        getPositionalReads(queryMat, algTibble()$count, algMetaData()$totalCount * 0.001, algMetaData()$maxCut1)
    })
    snp2Tibble <- reactive({
        queryMat <- query2Mat()
        refMat <- toupper(ref2Mat())
        queryMat[queryMat != refMat & queryMat != "-"] = "S"
        queryMat[queryMat == refMat] = "M"
        getPositionalReads(queryMat, algTibble()$count, algMetaData()$totalCount * 0.001, algMetaData()$maxCut2)
    })

    posBaseRef1TempFile <- tempfile(tmpdir=file.path("www", session$token))
    posBaseRef2TempFile <- tempfile(tmpdir=file.path("www", session$token))
    output$posBaseRef1Plot <- renderUI({
        req(input$algfiles)
        req(input$positionalMode)
        if (input$positionalMode == "histgram base") {
            drawPositionalStatic(base1Tibble(), insert1Count(), posBaseRef1TempFile)
        } else if (input$positionalMode == "histgram indel") {
            drawPositionalStatic(MSD1Tibble(), insert1Count(), posBaseRef1TempFile)
        } else if (input$positionalMode == "read base") {
            drawPositionalReads(read1Tibble(), algMetaData()$maxCut1, algMetaData()$maxCut1down, posBaseRef1TempFile)
        } else if (input$positionalMode == "read snp") {
            drawPositionalSnps(snp1Tibble(), algMetaData()$maxCut1, algMetaData()$maxCut1down, posBaseRef1TempFile)
        } else if (input$positionalMode == "logo probability") {
            drawPositionalLogo(base1Freq()[2:5,], "prob", "ACGT", posBaseRef1TempFile)
        } else if (input$positionalMode == "logo bits") {
            drawPositionalLogo(base1Freq()[2:5,], "bits", "ACGT", posBaseRef1TempFile)
        } else if (input$positionalMode == "logo custom") {
            drawPositionalLogo(base1Freq()[2:5,], "custom", "ACGT", posBaseRef1TempFile)
        }
    })
    output$posBaseRef2Plot <- renderUI({
        req(input$algfiles)
        req(input$positionalMode)
        if (input$positionalMode == "histgram base") {
            drawPositionalStatic(base2Tibble(), insert2Count(), posBaseRef2TempFile)
        } else if (input$positionalMode == "histgram indel") {
            drawPositionalStatic(MSD2Tibble(), insert2Count(), posBaseRef2TempFile)
        } else if (input$positionalMode == "read base") {
            drawPositionalReads(read2Tibble(), algMetaData()$maxCut2, algMetaData()$maxCut2down, posBaseRef2TempFile)
        } else if (input$positionalMode == "read snp") {
            drawPositionalSnps(snp2Tibble(), algMetaData()$maxCut2, algMetaData()$maxCut2down, posBaseRef2TempFile)
        } else if (input$positionalMode == "logo probability") {
            drawPositionalLogo(base2Freq()[2:5,], "prob", "ACGT", posBaseRef2TempFile)
        } else if (input$positionalMode == "logo bits") {
            drawPositionalLogo(base2Freq()[2:5,], "bits", "ACGT", posBaseRef2TempFile)
        } else if (input$positionalMode == "logo custom") {
            drawPositionalLogo(base2Freq()[2:5,], "custom", "ACGT", posBaseRef2TempFile)
        }
    })

    #############################
    # micro homology
    #############################
    uniqTibble <- reactive({
        algTibble() |> select(refId, cut1, cut2) |> rownames_to_column(var = "index") |> mutate(index = as.integer(index)) |> filter(refId == as.integer(microRefId())) |> summarise(index = first(index), cut1 = first(cut1), cut2 = first(cut2))
    })
    mhTibble <- reactive({
        refSeq <- algTibble()$refNoGap[uniqTibble()$index]
        ref1Len <- algTibble()$ref1Len[uniqTibble()$index]
        getMicroHomologyTibble(
            substr(refSeq, 1, ref1Len),
            substr(refSeq, ref1Len + 1, nchar(refSeq)),
            uniqTibble()$cut1, uniqTibble()$cut2
        )
    })
    mhTibbleSub <- reactive({
        mhTibble() |> filter(pos1up - pos1low >= input$microThres)
    })
    refEnd1Start2Tibble <- reactive({
        getRefEnd1Start2Tibble(algTibble(), as.integer(microRefId()))    
    })
    refEnd1Start2TibbleMicro <- reactive({
        getRefEnd1Start2TibbleMicro(refEnd1Start2Tibble(), mhTibbleSub())
    })

    # Use microRefId as a proxy of input$microRefId. updateSelectInput will not update input$microRefId until the client send the updated values back to the server. microRefId helps to block renderPlot until input$microRefId got updated.
    microRefId <- reactiveVal()
    observe({
        microRefId(input$microRefId)
    })
    observe({
        microRefId(NULL)
        updateSelectInput(inputId = "microRefId", choices = algTibble()$refId |> unique())
    }) |> bindEvent(input$algfiles)

    mhMatrixTempFile <- tempfile(tmpdir=file.path("www", session$token))
    output$mhMatrixPlot <- renderUI({
        req(input$algfiles)
        req(microRefId())
        cat("render plot")
        drawMicroHomologyHeatmap(mhTibbleSub(), refEnd1Start2TibbleMicro(), algMetaData()$maxCut1, algMetaData()$maxCut2, algMetaData()$maxCut1down, algMetaData()$maxCut2down, input$microMode, mhMatrixTempFile)
    })

    ##############################
    # classic classification
    ##############################
    indelTypeTibble <- reactive({
        getIndelTypes(algTibble())
    })
    indelTypeTibbleEx <- reactive({
        getIndelTypesEx(algTibble())
    })

    classifyTempFile <- tempfile(tmpdir=file.path("www", session$token))
    output$claClaPlot <- renderUI({
        req(input$algfiles)
        if (input$claClaDistinctTemp) {
            if (input$claClaMode == "pie") {
                indelTypeTibbleEx() |> indelTypePiePlot(classifyTempFile)
            } else if (input$claClaMode == "waffle") {
                indelTypeTibbleEx() |> indelTypeWafflePlot(classifyTempFile)
            }
        } else {
            if (input$claClaMode == "pie") {
                indelTypeTibble() |> indelTypePiePlot(classifyTempFile)
            } else if (input$claClaMode == "waffle") {
                indelTypeTibble() |> indelTypeWafflePlot(classifyTempFile)
            }
        }
    })

    ##############################
    # distribution plot
    ##############################
    distriTempFile <- tempfile(tmpdir=file.path("www", session$token))
    output$distriPlot <- renderUI({
        req(input$algfiles)
        req(input$distriTarget)
        if (input$distriMode == "discrete") {
            algTibble() |> mutate(randInsert = nchar(randInsert)) |> select(c("count", input$distriTarget)) |> discreteDistribution(distriTempFile)
        } else if (input$distriMode == "continuous") {
            algTibble() |> mutate(randInsert = nchar(randInsert)) |> select(c("count", input$distriTarget)) |> continuousDistribution(distriTempFile)
        }
    })

    ##############################
    # pairwise plot
    ##############################
    pairwiseTempFile <- tempfile(tmpdir=file.path("www", session$token), fileext=".pdf")
    output$pairwisePlot <- renderUI({
        req(input$algfiles)
        req(input$pairwiseX != input$pairwiseY)
        ggFig <- algTibble() |> mutate(randInsert = nchar(randInsert)) |> select(input$pairwiseX, input$pairwiseY) |> pairwisePlot(input$pairwiseXscale, input$pairwiseYscale, input$pairwiseMethod, input$pairwiseSpan)
        output$pairwiseWarning <- renderText({
            capture.output(ggFig, type = "message")
        })
        ggsave(pairwiseTempFile, plot = ggFig)
        tags$iframe(src = sub("^www/", "", pairwiseTempFile), height = "1200px", width = "100%")
    })

    ###############################
    # polygon insertion
    ###############################
    polyInsTibble <- reactive({
        getPolyInsTibble(algTibble())
    })
    polyInsTibble1 <- reactive({
        polyInsTibble() |> select(count, insPos1, insLen1) |> rename(insPos = insPos1, insLen = insLen1) |> unnest(c(insPos, insLen)) |> summarise(count = sum(count), .by = c(insPos, insLen))
    })
    polyInsTibble2 <- reactive({
        polyInsTibble() |> select(count, insPos2, insLen2) |> rename(insPos = insPos2, insLen = insLen2) |> unnest(c(insPos, insLen)) |> summarise(count = sum(count), .by = c(insPos, insLen))
    })
    polyXY1 <- reactive({
        getPolyXY(polyInsTibble1(), "down")
    })
    polyXY2 <- reactive({
        getPolyXY(polyInsTibble2(), "up")
    })

    polyInsert1TempFile <- tempfile(tmpdir=file.path("www", session$token))
    polyInsert2TempFile <- tempfile(tmpdir=file.path("www", session$token))
    output$polyInsert1Plot <- renderUI({
        req(input$algfiles)
        plotPolyInsTibble(polyXY1(), c(-algMetaData()$maxCut1, algMetaData()$maxCut1down + algMetaData()$maxRandInsert), polyInsert1TempFile)
    })
    output$polyInsert2Plot <- renderUI({
        req(input$algfiles)
        plotPolyInsTibble(polyXY2(), c(-algMetaData()$maxCut2 - algMetaData()$maxRandInsert, algMetaData()$maxCut2down), polyInsert2TempFile)
    })

    ###############################
    # arc deletion
    ###############################
    arcDelTibble <- reactive({
        getArcDelTibble(algTibble())
    })
    arcDelTibble1 <- reactive({
        arcDelTibble() |> select(count, delStart1, delEnd1) |> rename(delStart = delStart1, delEnd = delEnd1) |> unnest(c(delStart, delEnd)) |> summarise(count = sum(count), .by = c(delStart, delEnd))
    })
    arcDelTibble2 <- reactive({
        arcDelTibble() |> select(count, delStart2, delEnd2) |> rename(delStart = delStart2, delEnd = delEnd2) |> unnest(c(delStart, delEnd)) |> summarise(count = sum(count), .by = c(delStart, delEnd))
    })

    arcDelete1TempFile <- tempfile(tmpdir=file.path("www", session$token))
    arcDelete2TempFile <- tempfile(tmpdir=file.path("www", session$token))
    output$arcDelete1Plot <- renderUI({
        req(input$algfiles)
        plotArcDelTibble(arcDelTibble1(), c(-algMetaData()$maxCut1, algMetaData()$maxCut1down), arcDelete1TempFile)
    })
    output$arcDelete2Plot <- renderUI({
        req(input$algfiles)
        plotArcDelTibble(arcDelTibble2(), c(-algMetaData()$maxCut2, algMetaData()$maxCut2down), arcDelete2TempFile)
    })

    ################################
    # kpLogo
    ################################
    # use kpLogoRegion as a proxy of input$kpLogoRegion to prevent render output$kpLogoPlot twice (one before input$kpLogoRegion is returned by client, and one after that)
    kpLogoRegion <- reactiveVal(c(0, 0))
    observe({
        kpLogoRegion(input$kpLogoRegion)
    })
    observe({
        kpLogoRegion(c(0, 0))
        sgLen <- nchar(sgRNAs()[1])
        updateSliderInput(inputId = "kpLogoRegion", value = c(max(sgLen - 5, 1), sgLen), min = 1, max = sgLen, step = 1)
    }) |> bindEvent(input$sgRNAfile$datapath)
    algTarget <- reactive({
        getKpLogoAlgTarget(algTibble(), sgRNAs(), editTarget(), input$kpLogoCountThres, input$kpLogoMethod)
    })

    outputKpLogoTempFile <- tempfile(tmpdir=file.path("www", session$token))
    weightKpLogoTempFile <- tempfile(tmpdir=file.path("www", session$token))
    targetKpLogoTempFile <- tempfile(tmpdir=file.path("www", session$token))
    bgFileKpLogoTempFile <- tempfile(tmpdir=file.path("www", session$token))
    output$kpLogoIframe <- renderUI({
        req(input$algfiles)
        req(input$sgRNAfile)
        req(kpLogoRegion()[2] > 0)
        plotKpLogoAlgTarget(algTarget(), input$kpLogoMethod, kpLogoRegion(), input$kpLogoKmer, outputKpLogoTempFile, weightKpLogoTempFile, targetKpLogoTempFile, bgFileKpLogoTempFile)
    })

    #####################################
    # kmer frequencies
    #####################################
    # use kmerRange as a proxy of input$kmerRange to prevent render output$kmerIframe twice (one before input$kmerRange is returned by client, and one after that)
    kmerRange <- reactiveVal(c(0, 0))
    observe({
        kmerRange(input$kmerRange)
    })
    output$kmerRangeUI <- renderUI({
        req(input$sgRNAfile)
        kmerRange(c(0, 0))
        sgLen <- nchar(sgRNAs()[1])
        numericRangeInput("kmerRange", "kmer range", value = rep(max(1, sgLen - 3), 2), min = 1, max = sgLen, step = 1)
    })
    kmerTibble <- reactive({
        algTibble() |> mutate(kmer = substr(sgRNAs()[refId], kmerRange()[1], kmerRange()[2]), target = editTarget()) |> summarise(count = sum(count), .by = c(kmer, target))
    })

    kmerPdfTempFile <- tempfile(tmpdir=file.path("www", session$token))
    output$kmerIframe <- renderUI({
        req(input$algfiles)
        req(input$sgRNAfile)
        req(kmerRange()[2] > 0)
        plotKmerFrequencies(kmerTibble(), kmerPdfTempFile)
    })
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

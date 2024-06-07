library(shiny)
library(bslib)
library(tidyverse)
library(shinyWidgets)
library(reshape2)
library(scales)
library(ggseqlogo)
library(waffle)

options(shiny.maxRequestSize = 100 * 1024^3)

source("helpers/alignBrowser.R")
source("helpers/alignOne.R")
source("helpers/baseSubstitute.R")
source("helpers/positionalStatistics.R")
source("helpers/microHomology.R")
source("helpers/classicClassify.R")

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
        selectInput("positionalMode", "display mode", choices = c("histgram base", "histgram indel", "read base", "read snp", "logo probability", "logo bits", "logo custom"), selected = "histgram base"),
        plotOutput("posBaseRef1Plot"),
        plotOutput("posBaseRef2Plot")
    ),
    tabPanel("micro homology",
        selectInput("microRefId", "reference id", choices = NULL),
        numericInput("microThres", "minimal micro homology length", value = 4, min = 1),
        plotOutput("mh_matrix_plot")
    ),
    tabPanel("classic classify",
        checkboxInput("claClaDistinctTemp", "discriminate templated insertion and random insertion"),
        selectInput("claClaMode", "mode", choices = c("pie", "waffle")),
        plotOutput("claClaPlot", height = "1000px")
    ),
    tabPanel("continuous distribution",
        selectInput("conDistriTarget", "what", choices = c("alignment score", "cut1", "cut2", "ref1 length", "ref2 Length", "alignment end in ref1", "alignment start in ref2", "upstream insertion", "downstream insertion", "random insertion", "templated insertion", "insertion", "upstream deletion", "downstream deletion", "deletion"))
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
    headers <- reactive({
        algLines()[seq(1, length(algLines()), 3)] |> strsplit("\t")
    })
    counts <- reactive({
        headers() |> vapply(function(x) as.integer(x[2]), 0)
    })
    totalCount <- reactive({
        sum(counts())
    })
    scores <- reactive({
        headers() |> vapply(function(x) as.double(x[3]), 0)
    })
    refIds <- reactive({
        headers() |> vapply(function(x) as.integer(x[4]), 0)
    })
    refAllSeqs <- reactive({
        algLines()[seq(2, length(algLines()), 3)] |> str_replace_all("-", "")
    })
    cuts1 <- reactive({
        headers() |> vapply(function(x) as.integer(x[16]), 0)
    })
    ref1Lens <- reactive({
        ref1Lens <- refAllSeqs() |> vapply(function(x) gregexpr("[acgtn]", x)[[1]][2], 0)
        attr(ref1Lens, "names") <- NULL
        return(ref1Lens)
    })
    cuts2 <- reactive({
        headers() |> vapply(function(x) as.integer(x[17]), 0) - ref1Lens()
    })
    ref2Lens <- reactive({
        nchar(refAllSeqs()) - ref1Lens()
    })
    ref1ends <- reactive({
        headers() |> vapply(function(x) as.integer(x[8]), 0)
    })
    ref2starts <- reactive({
        headers() |> vapply(function(x) as.integer(x[11]), 0) - ref1Lens()
    })
    maxCut1 <- reactive({
        max(cuts1())
    })
    maxCut2 <- reactive({
        max(cuts2())
    })
    maxCut1down <- reactive({
        max(ref1Lens() - cuts1())
    })
    maxCut2down <- reactive({
        max(ref2Lens() - cuts2())
    })
    upInserts <- reactive({
        ifelse(ref1ends() > cuts1(), ref1ends() - cuts1(), 0)
    })
    upDeletes <- reactive({
        ifelse(ref1ends() < cuts1(), cuts1() - ref1ends(), 0)
    })
    randomInserts <- reactive({
        headers() |> vapply(function(x) as.integer(nchar(x[10])), 0)
    })
    downInserts <- reactive({
        ifelse(ref2starts() < cuts2(), cuts2() - ref2starts(), 0)
    })
    downDeletes <- reactive({
        ifelse(ref2starts() > cuts2(), ref2starts() - cuts2(), 0)
    })
    templatedInserts <- reactive({
        upInserts() + downInserts()
    })
    deletes <- reactive({
        upDeletes() + downDeletes()
    })
    inserts <- reactive({
        randomInserts() + templatedInserts()
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
    refGapMat <- reactive({
        algLines()[seq(2, length(algLines()), 3)] |> extendToSameLength()
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
    base1Freq <- reactive({
        getPositionalBaseFreq(query1Mat(), counts())
    })
    base2Freq <- reactive({
        getPositionalBaseFreq(query2Mat(), counts())
    })
    base1Tibble <- reactive({
        base1Freq() |> posMatrixToTibble(maxCut1())
    })
    base2Tibble <- reactive({
        base2Freq() |> posMatrixToTibble(maxCut2())
    })
    MSD1Tibble <- reactive({
        getPositionalMSDTibble(query1Mat(), ref1Mat(), counts(), maxCut1())
    })
    MSD2Tibble <- reactive({
        getPositionalMSDTibble(query2Mat(), ref2Mat(), counts(), maxCut2())
    })
    insert1Count <- reactive({
        refList <- algLines()[seq(2, length(algLines()), 3)] |> vapply(function(x) substr(x, 1, gregexpr("[acgtn]", x)[[1]][3] - 1), "") |> strsplit("")
        calInsertionCount(refList, cuts1(), maxCut1down())
    })
    insert2Count <- reactive({
        refList <- algLines()[seq(2, length(algLines()), 3)] |> vapply(function(x) substr(x, gregexpr("[acgtn]", x)[[1]][2] + 1, nchar(x)), "") |> strsplit("")
        calInsertionCount(refList, cuts2(), maxCut2down())
    })
    read1Tibble <- reactive({
        getPositionalReads(query1Mat(), counts(), totalCount() * 0.001, maxCut1())
    })
    read2Tibble <- reactive({
        getPositionalReads(query2Mat(), counts(), totalCount() * 0.001, maxCut2())
    })
    snp1Tibble <- reactive({
        queryMat <- query1Mat()
        refMat <- toupper(ref1Mat())
        queryMat[queryMat != refMat & queryMat != "-"] = "S"
        queryMat[queryMat == refMat] = "M"
        getPositionalReads(queryMat, counts(), totalCount() * 0.001, maxCut1())
    })
    snp2Tibble <- reactive({
        queryMat <- query2Mat()
        refMat <- toupper(ref2Mat())
        queryMat[queryMat != refMat & queryMat != "-"] = "S"
        queryMat[queryMat == refMat] = "M"
        getPositionalReads(queryMat, counts(), totalCount() * 0.001, maxCut2())
    })

    output$posBaseRef1Plot <- renderPlot({
        req(input$algfiles)
        req(input$positionalMode)

        if (input$positionalMode == "histgram base") {
            drawPositionalStatic(base1Tibble(), insert1Count())
        }
        else if (input$positionalMode == "histgram indel") {
            drawPositionalStatic(MSD1Tibble(), insert1Count())
        } else if (input$positionalMode == "read base") {
            drawPositionalReads(read1Tibble(), maxCut1(), maxCut1down())
        } else if (input$positionalMode == "read snp") {
            drawPositionalSnps(snp1Tibble(), maxCut1(), maxCut1down())
        } else if (input$positionalMode == "logo probability") {
            drawPositionalLogo(base1Freq()[2:5,], 'prob', "ACGT")
        } else if (input$positionalMode == "logo bits") {
            drawPositionalLogo(base1Freq()[2:5,], 'bits', "ACGT")
        } else if (input$positionalMode == "logo custom") {
            drawPositionalLogo(base1Freq()[2:5,], 'custom', "ACGT")
        }
    })

    output$posBaseRef2Plot <- renderPlot({
        req(input$algfiles)
        req(input$positionalMode)

        if (input$positionalMode == "histgram base") {
            drawPositionalStatic(base2Tibble(), insert2Count())
        }
        else if (input$positionalMode == "histgram indel") {
            drawPositionalStatic(MSD2Tibble(), insert2Count())
        } else if (input$positionalMode == "read base") {
            drawPositionalReads(read2Tibble(), maxCut2(), maxCut2down())
        } else if (input$positionalMode == "read snp") {
            drawPositionalSnps(snp2Tibble(), maxCut2(), maxCut2down())
        } else if (input$positionalMode == "logo probability") {
            drawPositionalLogo(base2Freq()[2:5,], 'prob', "ACGT")
        } else if (input$positionalMode == "logo bits") {
            drawPositionalLogo(base2Freq()[2:5,], 'bits', "ACGT")
        } else if (input$positionalMode == "logo custom") {
            drawPositionalLogo(base2Freq()[2:5,], 'custom', "ACGT")
        }
    })

    #############################
    # micro homology
    #############################
    uniqTibble <- reactive({
        tibble(refId = refIds(), index = seq_len(length(refIds())), cut1 = cuts1(), cut2 = cuts2()) |> filter(refId == as.integer(input$microRefId)) |> summarise(index = first(index), cut1 = first(cut1), cut2 = first(cut2))
    })
    observe({
        updateSelectInput(inputId = "microRefId", label = "reference id", choices = refIds() |> unique())
    }) |> bindEvent(input$algfiles)
    mhTibble <- reactive({
        refSeq <- refAllSeqs()[uniqTibble()$index]
        ref1Len <- ref1Lens()[uniqTibble()$index]
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
        getRefEnd1Start2Tibble(ref1ends(), ref2starts(), cuts1(), cuts2(), refIds(), counts(), as.integer(input$microRefId))    
    })
    refEnd1Start2TibbleMicro <- reactive({
        getRefEnd1Start2TibbleMicro(refEnd1Start2Tibble(), mhTibbleSub())
    })
    output$mh_matrix_plot <- renderPlot({
        req(input$algfiles)
        req(input$microRefId)
        drawMicroHomologyHeatmap(mhTibbleSub(), refEnd1Start2TibbleMicro(), maxCut1(), maxCut2(), maxCut1down(), maxCut2down())
    })

    ##############################
    # classic classification
    ##############################
    indelTypeTibble <- reactive({
        getIndelTypes(inserts(), deletes(), counts())
    })
    indelTypeTibbleEx <- reactive({
        getIndelTypesEx(templatedInserts(), deletes(), randomInserts(), counts())
    })

    output$claClaPlot <- renderPlot({
        req(input$algfiles)
        if (input$claClaDistinctTemp) {
            if (input$claClaMode == "pie") {
                indelTypePiePlot(indelTypeTibbleEx())
            } else if (input$claClaMode == "waffle") {
                indelTypeWafflePlot(indelTypeTibbleEx())
            }
        } else {
            if (input$claClaMode == "pie") {
                indelTypePiePlot(indelTypeTibble())
            } else if (input$claClaMode == "waffle") {
                indelTypeWafflePlot(indelTypeTibble())
            }
        }
    })
}

# Run the app ----
shinyApp(ui = ui, server = server)
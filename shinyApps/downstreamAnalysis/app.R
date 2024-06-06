library(shiny)
library(bslib)
library(tidyverse)
library(shinyWidgets)
library(reshape2)
library(scales)

options(shiny.maxRequestSize = 100 * 1024^3)

source("helpers/alignBrowser.R")
source("helpers/alignOne.R")
source("helpers/baseSubstitute.R")
source("helpers/positionalStatistics.R")
source("helpers/microHomology.R")

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
        headers() |> vapply(function(x) as.integer(x[11]), 0)
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

    output$posBaseRef1Plot <- renderPlot({
        req(input$algfiles)
        req(input$positionalMode)

        if (input$positionalMode == "histgram base") {
            drawPositionalStatic(base1Tibble(), insert1Count())
        }
        else if (input$positionalMode == "histgram indel") {
            drawPositionalStatic(MSD1Tibble(), insert1Count())
        } else if (input$positionalMode == "read base") {
            drawPositionalReads(query1Mat(), maxCut1())
        } else if (input$positionalMode == "read snp") {
            drawPositionalSnps(query1Mat(), ref1Mat(), maxCut1())
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
            drawPositionalReads(query2Mat(), maxCut2())
        } else if (input$positionalMode == "read snp") {
            drawPositionalSnps(query2Mat(), ref2Mat(), maxCut2())
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
        ref1 <- substr(refSeq, 1, ref1Len)
        ref2 <- substr(refSeq, ref1Len + 1, nchar(refSeq))
        getMicroHomologyTibble(ref1, ref2, uniqTibble()$cut1, uniqTibble()$cut2)
    })
    mhTibbleSub <- reactive({
        mhTibble() |> filter(pos1up - pos1low >= input$microThres)
    })
    refEnd1Start2Tibble <- reactive({
        tibble(
            pos1 = ref1ends() - cuts1(),
            pos2 = ref2starts() - ref1Lens() - cuts2(),
            refId = refIds(),
            count = counts()
        ) |> filter(refId == as.integer(input$microRefId)) |> summarise(count = sum(count), .by = c("pos1", "pos2")) |> mutate(shift = pos1 - pos2)
    })
    refEnd1Start2TibbleMicro <- reactive({
        joinTibble <- refEnd1Start2Tibble() |> left_join(mhTibbleSub(), by = "shift", relationship = "many-to-many")
        outRangeTibble <- joinTibble |> summarise(inRange = any(pos1 >= pos1low & pos1 <= pos1up), count = first(count), .by = c("pos1", "pos2")) |> filter(is.na(inRange) | !inRange) |> mutate(inRange = NULL, cls = 0)
        inRangeTibble <- joinTibble |> filter(pos1 >= pos1low, pos1 <= pos1up)
        if (nrow(inRangeTibble) == 0) {
            return(outRangeTibble)
        }
        inRangeTibble |> rowwise() |> mutate(pos1 = list(seq(pos1low, pos1up))) |> ungroup() |> unnest(pos1) |> mutate(pos2 = pos1 - shift) |> select(pos1, pos2, count, cls) |> bind_rows(outRangeTibble)
    })
    output$mh_matrix_plot <- renderPlot({
        req(input$algfiles)
        req(input$microRefId)
        mhPosTibble <- mhTibbleSub() |> mutate(nacol = NA) |> pivot_longer(cols = c("pos1low", "pos1up", "nacol"), values_to = "mhPos1") |> mutate(mhPos2 = mhPos1 - shift) |> select(mhPos1, mhPos2)
        log10p1 = trans_new("log10p1", function(x) log10(x + 1), function(x) 10^x - 1)
        ggplot(refEnd1Start2TibbleMicro()) +
            geom_tile(aes(x = pos2, y = pos1, fill = count), height = 1, width = 1) +
            geom_path(aes(x = mhPos2, y = mhPos1), data = mhPosTibble) +
            scale_x_continuous(limits = c(-maxCut2() - 1, maxCut2down() + 1), expand=c(0, 0)) +
            scale_y_continuous(limits = c(-maxCut1() - 1, maxCut1down() + 1), expand=c(0, 0)) +
            scale_fill_gradientn(limits = c(0, NA), colors = c("white", "red"), trans = log10p1) +
            scale_size_area(max_size = 2) +
            coord_equal(ratio = 1)
    })
}

# Run the app ----
shinyApp(ui = ui, server = server)
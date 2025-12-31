library(shiny)
library(bslib)
library(tidyverse)
library(rearrVis, lib.loc = "./library")
library(ggseqlogo, lib.loc = "./library")

options(shiny.maxRequestSize = 10 * 1024^3)

#####################
# frontend
#####################
ui <- page_sidebar(
  sidebar = sidebar(
    tooltip(
      fileInput("algfiles", "Alignment files", multiple = TRUE),
      ".alg or .alg.gz files, mutiple files are concatenated"
    ),
  ),
  navbarPage(
    title = NULL,
    tabPanel(
      title = tooltip(
        "browser",
        "a browser for alignment results"
      ),
      tooltip(
        numericInput(
          "browserDisplayCount",
          "Display count",
          value = 100,
          min = 1,
        ),
        "read count to display"
      ),
      tooltip(
        htmlOutput(
          "alignments",
          style = "font-family: Courier,courier; white-space: nowrap;",
          inline = TRUE
        ),
        "alignments"
      )
    ),
    tabPanel(
      title = tooltip(
        "base",
        "summarize frequency of base substitution"
      ),
      tooltip(
        uiOutput("baseSubFreqPlot"),
        "frequency of base substitution"
      )
    ),
    tabPanel(
      title = tooltip(
        "position",
        "postional statistics for alignments"
      ),
      tooltip(
        selectInput(
          "positionalMode",
          "display mode",
          choices = c(
            "histgram base",
            "histgram indel",
            "read base A",
            "read base C",
            "read base G",
            "read base T",
            "read match",
            "read snp",
            "logo probability",
            "logo bits",
            "logo custom"
          ),
          selected = "histgram base"
        ),
        "display mode of positional plot"
      ),
      tooltip(
        uiOutput("posBaseRef1Plot"),
        "positional plot of reference1"
      ),
      tooltip(
        uiOutput("posBaseRef2Plot"),
        "positional plot of reference2"
      )
    ),
    tabPanel(
      title = tooltip(
        "MMEJ",
        "display micro-homology diagram"
      ),
      tooltip(
        selectInput("microRefId", "reference id", choices = NULL),
        "reference index to analyze micro-homology"
      ),
      tooltip(
        selectInput(
          "microMode",
          "display mode",
          choices = c("repeat", "separate")
        ),
        "!repeat count over all microhomology positions, or !separate count evenly"
      ),
      tooltip(
        numericInput(
          "microThres",
          "micro-homology length threshold",
          value = 4,
          min = 1
        ),
        "minimal micro-homology length to display"
      ),
      tooltip(
        uiOutput("mhMatrixPlot"),
        "micro-homology diagram"
      )
    ),
    tabPanel(
      title = tooltip(
        "classify",
        r"(plot pie diagram of mutation types)"
      ),
      tooltip(
        checkboxInput(
          "claClaDistinctTemp",
          r"(discriminate templated\random insertion)"
        ),
        "discriminate templated\random insertion"
      ),
      tooltip(
        uiOutput("claClaPlot"),
        r"(pie diagram of mutation types)"
      )
    ),
    tabPanel(
      title = tooltip(
        "distribution",
        "plot histogram of alignment properties"
      ),
      tooltip(
        selectInput(
          "distriTarget",
          "what",
          choices = c(
            "score",
            "cut1",
            "cut2",
            "ref1Len",
            "ref2Len",
            "ref1End",
            "ref2Start",
            "upInsert",
            "downInsert",
            "randInsert",
            "templatedInsert",
            "insert",
            "upDelete",
            "downDelete",
            "delete"
          ),
          multiple = TRUE
        ),
        "properties to hist, multiple properties are supported"
      ),
      tooltip(
        selectInput(
          "distriMode",
          "mode",
          choices = c("continuous", "discrete")
        ),
        r"(use continuous\discrete histogram)"
      ),
      tooltip(
        uiOutput("distriPlot"),
        "histogram of alignment properties"
      )
    ),
    tabPanel(
      title = tooltip(
        "pairwise",
        "plot pairwise relationship of selected property"
      ),
      tooltip(
        selectInput(
          "pairwiseX",
          "x",
          choices = c(
            "score",
            "cut1",
            "cut2",
            "ref1Len",
            "ref2Len",
            "ref1End",
            "ref2Start",
            "upInsert",
            "downInsert",
            "randInsert",
            "templatedInsert",
            "insert",
            "upDelete",
            "downDelete",
            "delete"
          )
        ),
        "x-axis property"
      ),
      tooltip(
        selectInput(
          "pairwiseY",
          "y",
          choices = c(
            "score",
            "cut1",
            "cut2",
            "ref1Len",
            "ref2Len",
            "ref1End",
            "ref2Start",
            "upInsert",
            "downInsert",
            "randInsert",
            "templatedInsert",
            "insert",
            "upDelete",
            "downDelete",
            "delete"
          )
        ),
        "y-axis property"
      ),
      tooltip(
        selectInput(
          "pairwiseXscale",
          "x scale",
          choices = c("identity", "log10")
        ),
        r"(choose identity\log10 scale for x-axis)"
      ),
      tooltip(
        selectInput(
          "pairwiseYscale",
          "y scale",
          choices = c("identity", "log10")
        ),
        r"(choose identity\log10 scale for y-axis)"
      ),
      tooltip(
        selectInput(
          "pairwiseMethod",
          "method",
          choices = c("auto", "lm", "glm", "gam", "loess")
        ),
        "smooth method"
      ),
      tooltip(
        numericInput(
          "pairwiseSpan",
          "span",
          0.75,
          min = 0,
          max = 1,
          step = 0.01
        ),
        "span for smooth method"
      ),
      tooltip(
        uiOutput("pairwisePlot"),
        "pairwise relationship of selected property"
      )
    ),
    tabPanel(
      title = tooltip(
        "polygon",
        "plot insertion triangle diagram, height for count, width for insertion length"
      ),
      tooltip(
        uiOutput("polyInsert1Plot"),
        "insertion triangle diagram of reference1"
      ),
      tooltip(
        uiOutput("polyInsert2Plot"),
        "insertion triangle diagram of reference2"
      )
    ),
    tabPanel(
      title = tooltip(
        "arc",
        "deletion arc diagram, weight for count, width for deletion length"
      ),
      tooltip(
        uiOutput("arcDelete1Plot"),
        "deletion arc diagram of reference1"
      ),
      tooltip(
        uiOutput("arcDelete2Plot"),
        "deletion arc diagram of reference2"
      )
    ),
    tabPanel(
      title = tooltip(
        "kmer",
        "plot kmer frequency in sgRNA"
      ),
      tooltip(
        fileInput("sgRNAfile", "sgRNA file"),
        "sgRNA file"
      ),
      tooltip(
        selectInput(
          "editTarget",
          "edit target",
          choices = c(
            "templated insertion",
            "random insertion",
            "insertion",
            "deletion",
            "templated indel",
            "indel",
            "wild type"
          )
        ),
        "targeted sequences mutation type"
      ),
      tooltip(
        numericInput("kmerLowBound", "kmer lower bound", value = 18, min = 1),
        "lower bound range of sgRNA to count kmer"
      ),
      tooltip(
        numericInput("kmerUpBound", "kmer upper bound", value = 20, min = 1),
        "upper bound range of sgRNA to count kmer"
      ),
      tooltip(
        uiOutput("kmerIframe"),
        "kmer frequency in sgRNA"
      )
    )
  )
)

#################
# backend
#################
server <- function(input, output, session) {
  ################################
  # session start
  ################################
  file.path("www", session$token) |> dir.create(recursive = TRUE)
  ################################
  # session end
  ################################
  session$onSessionEnded(function() {
    unlink(file.path("www", session$token), recursive = TRUE)
  })
  ################################
  # sidebar
  ################################
  algTibble <- reactive({
    algLines <- lapply(
      input$algfiles$datapath,
      function(datapath) {
        readLines(gzfile(datapath))
      }
    ) |>
      unlist()
    algLines[seq(1, length(algLines), 3)] |>
      I() |>
      read_tsv(
        col_names = c(
          "index",
          "count",
          "score",
          "refId",
          "upDangle",
          "ref1Start",
          "query1Start",
          "ref1End",
          "query1End",
          "randInsert",
          "ref2Start",
          "query2Start",
          "ref2End",
          "query2End",
          "downDangle",
          "cut1",
          "cut2"
        ),
        col_types = "iidiciiiiciiiicii",
        col_select = -1,
        na = character()
      ) |>
      mutate(
        refLine = algLines[seq(2, length(algLines), 3)],
        queryLine = algLines[seq(3, length(algLines), 3)],
        refNoGap = refLine |> str_replace_all("-", ""),
        ref1Len = refNoGap |>
          vapply(
            function(x) gregexpr("[acgtn]", x)[[1]][2],
            0,
            USE.NAMES = FALSE
          ),
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

  # Use proxy$xxxxx as a proxy of input$xxxxx. update??????Input will not update input$xxxxx until the client send the updated values back to the server. proxy$xxxxx helps to block renderUI until input$xxxxx got updated.
  proxy <- reactiveValues()
  observe({
    for (name in names(proxy)) {
      proxy[[name]] <- input[[name]]
    }
  })

  ################################
  # alignment browser
  ################################
  output$alignments <- renderText({
    req(input$algfiles)
    req(input$browserDisplayCount)
    paste(
      rearrVis::getMarkdownFromAlign(algTibble()[
        seq_len(min(input$browserDisplayCount, nrow(algTibble()))),
      ]),
      collapse = ""
    )
  })

  #########################################################
  # base substitution frequency
  #########################################################
  baseSubFreqTempFile <- tempfile(
    tmpdir = file.path("www", session$token),
    fileext = ".pdf"
  )
  output$baseSubFreqPlot <- renderUI({
    req(input$algfiles)
    ggFig <- algTibble() |>
      rearrVis::countBaseSubstitute() |>
      ggplot(aes(sub, count)) +
      geom_col() +
      scale_y_continuous(expand = c(0, 0))
    ggsave(
      baseSubFreqTempFile,
      plot = ggFig,
      height = 1800,
      width = 3600,
      unit = "px"
    )
    tags$iframe(
      src = sub("^www/", "", baseSubFreqTempFile),
      height = "1200px",
      width = "100%"
    )
  })

  #################################
  # positional statistics
  #################################
  refGapMat <- reactive({
    algTibble()$refLine |> rearrVis::extendToSameLength()
  })
  ref1Mat <- reactive({
    algTibble()$refNoGap |>
      substr(1, algTibble()$ref1Len) |>
      rearrVis::extendToAlignCut(algTibble()$cut1)
  })
  ref2Mat <- reactive({
    algTibble()$refNoGap |>
      substr(
        algTibble()$ref1Len + 1,
        algTibble()$ref1Len + algTibble()$ref2Len
      ) |>
      rearrVis::extendToAlignCut(algTibble()$cut2)
  })
  queryAllSeqs <- reactive({
    tmpMat <- algTibble()$queryLine |>
      rearrVis::extendToSameLength() |>
      t() |>
      as.vector()
    rearrVis::vectorToStringVector(
      tmpMat[(refGapMat() |> t() |> as.vector()) != "-"],
      algTibble()$ref1Len + algTibble()$ref2Len
    )
  })
  query1Mat <- reactive({
    queryAllSeqs() |>
      substr(1, algTibble()$ref1Len) |>
      rearrVis::extendToAlignCut(algTibble()$cut1)
  })
  query2Mat <- reactive({
    queryAllSeqs() |>
      substr(algTibble()$ref1Len + 1, nchar(queryAllSeqs())) |>
      rearrVis::extendToAlignCut(algTibble()$cut2)
  })
  base1Freq <- reactive({
    rearrVis::getPositionalBaseFreq(query1Mat(), algTibble()$count)
  })
  base2Freq <- reactive({
    rearrVis::getPositionalBaseFreq(query2Mat(), algTibble()$count)
  })
  base1Tibble <- reactive({
    base1Freq() |> rearrVis::posMatrixToTibble(algMetaData()$maxCut1)
  })
  base2Tibble <- reactive({
    base2Freq() |> rearrVis::posMatrixToTibble(algMetaData()$maxCut2)
  })
  MSD1Tibble <- reactive({
    rearrVis::getPositionalMSDTibble(
      query1Mat(),
      ref1Mat(),
      algTibble()$count,
      algMetaData()$maxCut1
    )
  })
  MSD2Tibble <- reactive({
    rearrVis::getPositionalMSDTibble(
      query2Mat(),
      ref2Mat(),
      algTibble()$count,
      algMetaData()$maxCut2
    )
  })
  insert1Count <- reactive({
    algTibble()$refLine |>
      vapply(
        function(x) substr(x, 1, gregexpr("[acgtn]", x)[[1]][3] - 1),
        "",
        USE.NAMES = FALSE
      ) |>
      strsplit("") |>
      rearrVis::calInsertionCount(
        cuts = algTibble()$cut1,
        maxCutDown = algMetaData()$maxCut1down
      )
  })
  insert2Count <- reactive({
    algTibble()$refLine |>
      vapply(
        function(x) substr(x, gregexpr("[acgtn]", x)[[1]][2] + 1, nchar(x)),
        "",
        USE.NAMES = FALSE
      ) |>
      strsplit("") |>
      rearrVis::calInsertionCount(
        cuts = algTibble()$cut2,
        maxCutDown = algMetaData()$maxCut2down
      )
  })
  read1ATibble <- reactive({
    rearrVis::getPositionalReads(
      query1Mat() == "A",
      algTibble()$count,
      200,
      algMetaData()$maxCut1
    )
  })
  read2ATibble <- reactive({
    rearrVis::getPositionalReads(
      query2Mat() == "A",
      algTibble()$count,
      200,
      algMetaData()$maxCut2
    )
  })
  read1CTibble <- reactive({
    rearrVis::getPositionalReads(
      query1Mat() == "C",
      algTibble()$count,
      200,
      algMetaData()$maxCut1
    )
  })
  read2CTibble <- reactive({
    rearrVis::getPositionalReads(
      query2Mat() == "C",
      algTibble()$count,
      200,
      algMetaData()$maxCut2
    )
  })
  read1GTibble <- reactive({
    rearrVis::getPositionalReads(
      query1Mat() == "G",
      algTibble()$count,
      200,
      algMetaData()$maxCut1
    )
  })
  read2GTibble <- reactive({
    rearrVis::getPositionalReads(
      query2Mat() == "G",
      algTibble()$count,
      200,
      algMetaData()$maxCut2
    )
  })
  read1TTibble <- reactive({
    rearrVis::getPositionalReads(
      query1Mat() == "T",
      algTibble()$count,
      200,
      algMetaData()$maxCut1
    )
  })
  read2TTibble <- reactive({
    rearrVis::getPositionalReads(
      query2Mat() == "T",
      algTibble()$count,
      200,
      algMetaData()$maxCut2
    )
  })
  match1Tibble <- reactive({
    rearrVis::getPositionalReads(
      query1Mat() == toupper(ref1Mat()),
      algTibble()$count,
      200,
      algMetaData()$maxCut1
    )
  })
  match2Tibble <- reactive({
    rearrVis::getPositionalReads(
      query2Mat() == toupper(ref2Mat()),
      algTibble()$count,
      200,
      algMetaData()$maxCut2
    )
  })
  snp1Tibble <- reactive({
    rearrVis::getPositionalReads(
      query1Mat() != toupper(ref1Mat()) & query1Mat() != "-",
      algTibble()$count,
      200,
      algMetaData()$maxCut1
    )
  })
  snp2Tibble <- reactive({
    rearrVis::getPositionalReads(
      query2Mat() != toupper(ref2Mat()) & query2Mat() != "-",
      algTibble()$count,
      200,
      algMetaData()$maxCut2
    )
  })

  ref1TempFile <- tempfile(
    tmpdir = file.path("www", session$token),
    fileext = ".pdf"
  )
  ref2TempFile <- tempfile(
    tmpdir = file.path("www", session$token),
    fileext = ".pdf"
  )
  output$posBaseRef1Plot <- renderUI({
    req(input$algfiles)
    req(input$positionalMode)
    if (input$positionalMode == "histgram base") {
      rearrVis::drawPositionalStatic(
        base1Tibble(),
        insert1Count(),
        ref1TempFile
      )
    } else if (input$positionalMode == "histgram indel") {
      rearrVis::drawPositionalStatic(
        MSD1Tibble(),
        insert1Count(),
        ref1TempFile
      )
    } else if (input$positionalMode == "read base A") {
      rearrVis::drawPositionalReads(
        read1ATibble(),
        algMetaData()$maxCut1,
        algMetaData()$maxCut1down,
        ref1TempFile
      )
    } else if (input$positionalMode == "read base C") {
      rearrVis::drawPositionalReads(
        read1CTibble(),
        algMetaData()$maxCut1,
        algMetaData()$maxCut1down,
        ref1TempFile
      )
    } else if (input$positionalMode == "read base G") {
      rearrVis::drawPositionalReads(
        read1GTibble(),
        algMetaData()$maxCut1,
        algMetaData()$maxCut1down,
        ref1TempFile
      )
    } else if (input$positionalMode == "read base T") {
      rearrVis::drawPositionalReads(
        read1TTibble(),
        algMetaData()$maxCut1,
        algMetaData()$maxCut1down,
        ref1TempFile
      )
    } else if (input$positionalMode == "read match") {
      rearrVis::drawPositionalReads(
        match1Tibble(),
        algMetaData()$maxCut1,
        algMetaData()$maxCut1down,
        ref1TempFile
      )
    } else if (input$positionalMode == "read snp") {
      rearrVis::drawPositionalReads(
        snp1Tibble(),
        algMetaData()$maxCut1,
        algMetaData()$maxCut1down,
        ref1TempFile
      )
    } else if (input$positionalMode == "logo probability") {
      rearrVis::drawPositionalLogo(
        base1Freq()[2:5, ],
        "prob",
        "ACGT",
        ref1TempFile
      )
    } else if (input$positionalMode == "logo bits") {
      rearrVis::drawPositionalLogo(
        base1Freq()[2:5, ],
        "bits",
        "ACGT",
        ref1TempFile
      )
    } else if (input$positionalMode == "logo custom") {
      rearrVis::drawPositionalLogo(
        base1Freq()[2:5, ],
        "custom",
        "ACGT",
        ref1TempFile
      )
    }
  })
  output$posBaseRef2Plot <- renderUI({
    req(input$algfiles)
    req(input$positionalMode)
    if (input$positionalMode == "histgram base") {
      rearrVis::drawPositionalStatic(
        base2Tibble(),
        insert2Count(),
        ref2TempFile
      )
    } else if (input$positionalMode == "histgram indel") {
      rearrVis::drawPositionalStatic(
        MSD2Tibble(),
        insert2Count(),
        ref2TempFile
      )
    } else if (input$positionalMode == "read base A") {
      rearrVis::drawPositionalReads(
        read2ATibble(),
        algMetaData()$maxCut2,
        algMetaData()$maxCut2down,
        ref2TempFile
      )
    } else if (input$positionalMode == "read base C") {
      rearrVis::drawPositionalReads(
        read2CTibble(),
        algMetaData()$maxCut2,
        algMetaData()$maxCut2down,
        ref2TempFile
      )
    } else if (input$positionalMode == "read base G") {
      rearrVis::drawPositionalReads(
        read2GTibble(),
        algMetaData()$maxCut2,
        algMetaData()$maxCut2down,
        ref2TempFile
      )
    } else if (input$positionalMode == "read base T") {
      rearrVis::drawPositionalReads(
        read2TTibble(),
        algMetaData()$maxCut2,
        algMetaData()$maxCut2down,
        ref2TempFile
      )
    } else if (input$positionalMode == "read match") {
      rearrVis::drawPositionalReads(
        match2Tibble(),
        algMetaData()$maxCut2,
        algMetaData()$maxCut2down,
        ref2TempFile
      )
    } else if (input$positionalMode == "read snp") {
      rearrVis::drawPositionalReads(
        snp2Tibble(),
        algMetaData()$maxCut2,
        algMetaData()$maxCut2down,
        ref2TempFile
      )
    } else if (input$positionalMode == "logo probability") {
      rearrVis::drawPositionalLogo(
        base2Freq()[2:5, ],
        "prob",
        "ACGT",
        ref2TempFile
      )
    } else if (input$positionalMode == "logo bits") {
      rearrVis::drawPositionalLogo(
        base2Freq()[2:5, ],
        "bits",
        "ACGT",
        ref2TempFile
      )
    } else if (input$positionalMode == "logo custom") {
      rearrVis::drawPositionalLogo(
        base2Freq()[2:5, ],
        "custom",
        "ACGT",
        ref2TempFile
      )
    }
  })

  #############################
  # micro homology
  #############################
  uniqTibble <- reactive({
    algTibble() |>
      select(refId, cut1, cut2) |>
      rownames_to_column(var = "index") |>
      mutate(index = as.integer(index)) |>
      filter(refId == as.integer(proxy$microRefId)) |>
      summarise(index = first(index), cut1 = first(cut1), cut2 = first(cut2))
  })
  mhTibble <- reactive({
    refSeq <- algTibble()$refNoGap[uniqTibble()$index]
    ref1Len <- algTibble()$ref1Len[uniqTibble()$index]
    rearrVis::getMicroHomologyTibble(
      substr(refSeq, 1, ref1Len),
      substr(refSeq, ref1Len + 1, nchar(refSeq)),
      uniqTibble()$cut1,
      uniqTibble()$cut2
    )
  })
  mhTibbleSub <- reactive({
    mhTibble() |> filter(pos1up - pos1low >= input$microThres)
  })
  refEnd1Start2Tibble <- reactive({
    rearrVis::getRefEnd1Start2Tibble(algTibble(), as.integer(proxy$microRefId))
  })
  refEnd1Start2TibbleMicro <- reactive({
    rearrVis::getRefEnd1Start2TibbleMicro(refEnd1Start2Tibble(), mhTibbleSub())
  })

  observe({
    proxy$microRefId <- NULL
    updateSelectInput(
      inputId = "microRefId",
      choices = algTibble()$refId |> unique()
    )
  }) |>
    bindEvent(input$algfiles)

  mhMatrixTempFile <- tempfile(
    tmpdir = file.path("www", session$token),
    fileext = ".pdf"
  )
  output$mhMatrixPlot <- renderUI({
    req(input$algfiles)
    req(proxy$microRefId)
    rearrVis::drawMicroHomologyHeatmap(
      mhTibbleSub(),
      refEnd1Start2TibbleMicro(),
      algMetaData()$maxCut1,
      algMetaData()$maxCut2,
      algMetaData()$maxCut1down,
      algMetaData()$maxCut2down,
      input$microMode,
      mhMatrixTempFile
    )
  })

  ##############################
  # classic classification
  ##############################
  indelTypeTibble <- reactive({
    rearrVis::getIndelTypes(algTibble())
  })
  indelTypeTibbleEx <- reactive({
    rearrVis::getIndelTypesEx(algTibble())
  })

  classifyTempFile <- tempfile(
    tmpdir = file.path("www", session$token),
    fileext = ".pdf"
  )
  output$claClaPlot <- renderUI({
    req(input$algfiles)
    if (input$claClaDistinctTemp) {
      indelTypeTibbleEx() |> rearrVis::indelTypePiePlot(classifyTempFile)
    } else {
      indelTypeTibble() |> rearrVis::indelTypePiePlot(classifyTempFile)
    }
  })

  ##############################
  # distribution plot
  ##############################
  distriTempFile <- tempfile(
    tmpdir = file.path("www", session$token),
    fileext = ".pdf"
  )
  output$distriPlot <- renderUI({
    req(input$algfiles)
    req(input$distriTarget)
    if (input$distriMode == "discrete") {
      algTibble() |>
        mutate(randInsert = nchar(randInsert)) |>
        select(c("count", input$distriTarget)) |>
        rearrVis::discreteDistribution(distriTempFile)
    } else if (input$distriMode == "continuous") {
      algTibble() |>
        mutate(randInsert = nchar(randInsert)) |>
        select(c("count", input$distriTarget)) |>
        rearrVis::continuousDistribution(distriTempFile)
    }
  })

  ##############################
  # pairwise plot
  ##############################
  pairwiseTempFile <- tempfile(
    tmpdir = file.path("www", session$token),
    fileext = ".pdf"
  )
  output$pairwisePlot <- renderUI({
    req(input$algfiles)
    req(input$pairwiseX != input$pairwiseY)
    algTibble() |>
      mutate(randInsert = nchar(randInsert)) |>
      select(input$pairwiseX, input$pairwiseY) |>
      rearrVis::pairwisePlot(
        input$pairwiseXscale,
        input$pairwiseYscale,
        input$pairwiseMethod,
        input$pairwiseSpan,
        pairwiseTempFile
      )
  })

  ###############################
  # polygon insertion
  ###############################
  polyInsTibble <- reactive({
    rearrVis::getPolyInsTibble(algTibble())
  })
  polyInsTibble1 <- reactive({
    polyInsTibble() |>
      select(count, insPos1, insLen1) |>
      rename(insPos = insPos1, insLen = insLen1) |>
      unnest(c(insPos, insLen)) |>
      summarise(count = sum(count), .by = c(insPos, insLen))
  })
  polyInsTibble2 <- reactive({
    polyInsTibble() |>
      select(count, insPos2, insLen2) |>
      rename(insPos = insPos2, insLen = insLen2) |>
      unnest(c(insPos, insLen)) |>
      summarise(count = sum(count), .by = c(insPos, insLen))
  })
  polyXY1 <- reactive({
    rearrVis::getPolyXY(polyInsTibble1(), "down")
  })
  polyXY2 <- reactive({
    rearrVis::getPolyXY(polyInsTibble2(), "up")
  })

  polyInsert1TempFile <- tempfile(
    tmpdir = file.path("www", session$token),
    fileext = ".pdf"
  )
  polyInsert2TempFile <- tempfile(
    tmpdir = file.path("www", session$token),
    fileext = ".pdf"
  )
  output$polyInsert1Plot <- renderUI({
    req(input$algfiles)
    rearrVis::plotPolyInsTibble(
      polyXY1(),
      c(
        -algMetaData()$maxCut1,
        algMetaData()$maxCut1down + algMetaData()$maxRandInsert
      ),
      polyInsert1TempFile
    )
  })
  output$polyInsert2Plot <- renderUI({
    req(input$algfiles)
    rearrVis::plotPolyInsTibble(
      polyXY2(),
      c(
        -algMetaData()$maxCut2 - algMetaData()$maxRandInsert,
        algMetaData()$maxCut2down
      ),
      polyInsert2TempFile
    )
  })

  ###############################
  # arc deletion
  ###############################
  arcDelTibble <- reactive({
    rearrVis::getArcDelTibble(algTibble())
  })
  arcDelTibble1 <- reactive({
    arcDelTibble() |>
      select(count, delStart1, delEnd1) |>
      rename(delStart = delStart1, delEnd = delEnd1) |>
      unnest(c(delStart, delEnd)) |>
      summarise(count = sum(count), .by = c(delStart, delEnd))
  })
  arcDelTibble2 <- reactive({
    arcDelTibble() |>
      select(count, delStart2, delEnd2) |>
      rename(delStart = delStart2, delEnd = delEnd2) |>
      unnest(c(delStart, delEnd)) |>
      summarise(count = sum(count), .by = c(delStart, delEnd))
  })
  arcDelSegTibble1 <- reactive({
    rearrVis::getArcSegment(arcDelTibble1(), 100)
  })
  arcDelSegTibble2 <- reactive({
    rearrVis::getArcSegment(arcDelTibble2(), 100)
  })

  arcDelete1TempFile <- tempfile(
    tmpdir = file.path("www", session$token),
    fileext = ".pdf"
  )
  arcDelete2TempFile <- tempfile(
    tmpdir = file.path("www", session$token),
    fileext = ".pdf"
  )
  output$arcDelete1Plot <- renderUI({
    req(input$algfiles)
    rearrVis::plotArcDelTibble(
      arcDelSegTibble1(),
      c(-algMetaData()$maxCut1, algMetaData()$maxCut1down),
      arcDelete1TempFile
    )
  })
  output$arcDelete2Plot <- renderUI({
    req(input$algfiles)
    rearrVis::plotArcDelTibble(
      arcDelSegTibble2(),
      c(-algMetaData()$maxCut2, algMetaData()$maxCut2down),
      arcDelete2TempFile
    )
  })

  #####################################
  # kmer frequencies
  #####################################
  sgRNAs <- reactive({
    readLines(input$sgRNAfile$datapath)
  })

  editTarget <- reactive({
    if (input$editTarget == "templated insertion") {
      return(as.logical(algTibble()$templatedInsert))
    } else if (input$editTarget == "random insertion") {
      return(algTibble()$randInsert != "")
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

  kmerTibble <- reactive({
    req(input$kmerLowBound <= input$kmerUpBound)
    algTibble() |>
      mutate(
        kmer = substr(
          sgRNAs()[refId],
          input$kmerLowBound,
          min(input$kmerUpBound, sgRNAs() |> vapply(nchar, integer(1)) |> min())
        ),
        target = editTarget()
      ) |>
      summarise(count = sum(count), .by = c(kmer, target))
  })

  kmerPdfTempFile <- tempfile(
    tmpdir = file.path("www", session$token),
    fileext = ".pdf"
  )
  output$kmerIframe <- renderUI({
    req(input$algfiles)
    req(input$sgRNAfile)
    rearrVis::plotKmerFrequencies(kmerTibble(), kmerPdfTempFile)
  })
}

################################
# end app
################################
onStop(function() {
  unlink("www", recursive = TRUE)
})

################################
# run app
################################
shinyApp(ui = ui, server = server)

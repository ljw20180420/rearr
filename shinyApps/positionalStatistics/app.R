library(shiny)
library(bslib)
library(tidyverse)

options(shiny.maxRequestSize = 100 * 1024^2)

# Define UI ----
ui <- page_sidebar(
    title = "Positional Statistics",
    sidebar = sidebar(
        helpText("This APP show positional statistics of query reads relative to reference."),
        selectInput(
            "mode",
            "Display Mode",
            choices = c("base", "indel"),
            selected = "base"
        ),
        fileInput("algfile", label = "Alignment file"),
        numericInput("cut1", "cut1", value = NA, min = 0),
        numericInput("cut2", "cut2", value = NA, min = 0),
    ),
    plotOutput("posBaseRef1Plot"),
    plotOutput("posBaseRef2Plot")
)

# Define server logic ----
server <- function(input, output) {
    algLines <- reactive({
        req(input$algfile)
        readLines(input$algfile$datapath)
    })
    counts <- reactive({
        algLines()[seq(1, length(algLines()), 3)] |> strsplit("\t") |> vapply(function(x) as.integer(x[2]), 0)
    })
    ref <- reactive({
        str_replace_all(algLines()[2], "-", "")
    })
    ref1Len <- reactive({
        gregexpr("[acgtn]", ref())[[1]][2]
    })
    maxLen <- reactive({
        max(vapply(algLines()[seq(2, length(algLines()), 3)], nchar, 0, USE.NAMES=FALSE))
    })
    refMat <- reactive({
        algLines()[seq(2, length(algLines()), 3)] |> str_pad(width = maxLen(), side = 'right', pad = "-") |> strsplit("") |> unlist() |> matrix(ncol = maxLen(), byrow = TRUE)
    })
    queryMat <- reactive({
        req(input$algfile)
        tmpMat <- algLines()[seq(3, length(algLines()), 3)] |> str_pad(width = maxLen(), side = 'right', pad = "-") |> strsplit("") |> unlist() |> matrix(ncol = maxLen(), byrow = TRUE)
        tmpMat[refMat() != '-'] |> matrix(ncol = nchar(ref()))
    })
    baseFreq <- reactive({
        rbind(
            colSums((queryMat() == "A") * counts()),
            colSums((queryMat() == "C") * counts()),
            colSums((queryMat() == "G") * counts()),
            colSums((queryMat() == "T") * counts()),
            colSums((queryMat() == "-") * counts())
        )
    })
    MSDFreq <- reactive({
        tmpMat <- refMat()[refMat() != '-'] |> matrix(ncol = nchar(ref())) |> toupper()
        rbind(
            colSums((queryMat() == tmpMat) * counts()),
            colSums((queryMat() != tmpMat & queryMat() != '-') * counts()),
            colSums((queryMat() == '-') * counts())
        )
    })
    insertCount <- reactive({
        refList <- algLines()[seq(2, length(algLines()), 3)] |> strsplit("")
        insertList <- vector("list", length(refList))
        for (i in seq_len(length(refList))) {
            mask <- refList[[i]] != "-"
            insertList[[i]] <- cumsum(mask)[!mask]
        }
        histCount <- insertList |> unlist() |> table()
        insertCount <- rep(0, nchar(ref()) + 1)
        insertCount[as.integer(names(histCount)) + 1] <- histCount
        insertCount
    })

    output$posBaseRef1Plot <- renderPlot({
        req(input$mode)
        req(input$cut1)
        if (input$mode == "base")
        {
            tibble(count = c(baseFreq()), pos = rep(seq(ncol(baseFreq())) - 0.5, each = nrow(baseFreq())), base = rep(c("A", "C", "G", "T", "-"), times = ncol(baseFreq()))) |>
            filter(pos <= ref1Len()) |>
            mutate(base = factor(base, levels = c("-", "A", "C", "G", "T"))) |>
            mutate(rel1pos = pos - input$cut1) |>
            ggplot(aes(rel1pos, count)) +
            geom_col(aes(fill = base)) +
            geom_step(aes(pos - input$cut1, count, color = "black"), data = tibble(pos = 0:nchar(ref()), count = insertCount()) |> filter(pos > 0, pos < ref1Len()), direction = "mid") +
            scale_x_continuous(name = "position relative to cut1", expand = c(0,0)) +
            scale_y_continuous(expand = c(0,0)) +
            scale_color_identity(name = NULL, guide = guide_legend(), labels = "insertion")
        }
        else {
            tibble(count = c(MSDFreq()), pos = rep(seq(ncol(MSDFreq())) - 0.5, each = nrow(MSDFreq())), type = rep(c("match", "SNP", "delete"), times = ncol(MSDFreq()))) |>
                filter(pos <= ref1Len()) |>
                mutate(type = factor(type, levels = c("delete", "SNP", "match"))) |>
                mutate(rel1pos = pos - input$cut1) |>
                ggplot(aes(rel1pos, count)) +
                geom_col(aes(fill = type)) +
                geom_step(aes(pos - input$cut1, count), data = tibble(pos = 0:nchar(ref()), count = insertCount()) |> filter(pos > 0, pos < ref1Len()), color = "black", direction = "mid") +
                scale_x_continuous(name = "position relative to cut1", expand = c(0,0)) +
                scale_y_continuous(expand = c(0,0)) +
                scale_color_identity(name = NULL, guide = guide_legend(), labels = "insertion")
        }
    })

    output$posBaseRef2Plot <- renderPlot({
        req(input$mode)
        req(input$cut2)
        if (input$mode == "base") {
            tibble(count = c(baseFreq()), pos = rep(seq(ncol(baseFreq())) - 0.5, each = nrow(baseFreq())), base = rep(c("A", "C", "G", "T", "-"), times = ncol(baseFreq()))) |>
            filter(pos > ref1Len()) |>
            mutate(base = factor(base, levels = c("-", "A", "C", "G", "T"))) |>
            mutate(rel2pos = pos - ref1Len() - input$cut2) |>
            ggplot(aes(rel2pos, count)) +
            geom_col(aes(fill = base)) +
            geom_step(aes(pos - ref1Len() - input$cut2, count, color = "black"), data = tibble(pos = 0:nchar(ref()), count = insertCount()) |> filter(pos > ref1Len(), pos < nchar(ref())), direction = "mid") +
            scale_x_continuous(name = "position relative to cut2", expand = c(0,0)) +
            scale_y_continuous(expand = c(0,0)) +
            scale_color_identity(name = NULL, guide = guide_legend(), labels = "insertion")
        }
        else {
            tibble(count = c(MSDFreq()), pos = rep(seq(ncol(MSDFreq())) - 0.5, each = nrow(MSDFreq())), type = rep(c("match", "SNP", "delete"), times = ncol(MSDFreq()))) |>
                filter(pos > ref1Len()) |>
                mutate(type = factor(type, levels = c("delete", "SNP", "match"))) |>
                mutate(rel2pos = pos - ref1Len() - input$cut2) |>
                ggplot(aes(rel2pos, count)) +
                geom_col(aes(fill = type)) +
                geom_step(aes(pos - ref1Len() - input$cut2, count), data = tibble(pos = 0:nchar(ref()), count = insertCount()) |> filter(pos > ref1Len(), pos < nchar(ref())), color = "black", direction = "mid") +
                scale_x_continuous(name = "position relative to cut2", expand = c(0,0)) +
                scale_y_continuous(expand = c(0,0)) +
                scale_color_identity(name = NULL, guide = guide_legend(), labels = "insertion")
        }
    })
}

# Run the app ----
shinyApp(ui = ui, server = server)
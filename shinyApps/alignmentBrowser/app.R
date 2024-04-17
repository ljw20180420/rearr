library(shiny)
library(bslib)
library(tidyverse)

options(shiny.maxRequestSize = 100 * 1024^2)

arrangeInsertion <- function(query, insertion, insertionCollapse) {
    insertionLines <- c("")
    for (i in seq_len(length(insertionCollapse))) {
        insertionAdded <- FALSE
        insertionSeg <- substr(query, insertion[i], insertion[i] + attr(insertion, which = "match.length")[i] - 1)
        for (j in seq_len(length(insertionLines))) {
            if (nchar(insertionLines[j]) < insertionCollapse[i]) {
                insertionLines[j] <- sprintf("%s%s%s.",
                    insertionLines[j],
                    strrep(".", insertionCollapse[i] - nchar(insertionLines[j]) - 1),
                    insertionSeg
                )
                insertionAdded <- TRUE
                break
            }
        }
        if (!insertionAdded)
        {
            insertionLines <- c(
                insertionLines,
                sprintf("%s%s",
                    strrep(".", insertionCollapse[i] - 1),
                    insertionSeg
                ) 
            )
        }
    }
    # add a space at the beginning to fit the start bar trick of query
    insertionLines <- paste0(" ", insertionLines)
    return(insertionLines)
}

# Define UI ----
ui <- page_sidebar(
    tags$head(
        tags$style(HTML("
            .sequence {
                font-family: Courier,courier;
                white-space: nowrap;
                position: sticky;
                top: -24px;
                margin: 0, 0, 0, 0;
                padding: 0, 0, 0, 0;
                border: none;
            }
            .query {
                font-family: Courier,courier;
                white-space: nowrap;
            }
        "))
    ),

    # title = "Alignment Browser",
    sidebar = sidebar(
        helpText("This APP browser an alignment file."),
        fileInput("algfile", label = "Alignment file"),
        numericInput("cut1", "cut1", value = NA, min = 0),
        numericInput("cut2", "cut2", value = NA, min = 0),
    ),
    htmlOutput("reference", class = "sequence", inline = TRUE),
    htmlOutput("queries", class = "query", inline = TRUE)
)

# Define server logic ----
server <- function(input, output) {
    reference <- reactive({
        req(input$algfile)
        readLines(input$algfile$datapath, n = 2)[2] |> str_replace_all("[- ]", "")
    })
    ref1Len <- reactive({
        gregexpr("[acgtn]", reference())[[1]][2]
    })

    algLines <- reactive({
        req(input$algfile)
        readLines(input$algfile$datapath)
    })

    counts <- reactive({
        algLines()[seq(1, length(algLines()), 3)] |> strsplit("\t") |> vapply(function(x) as.integer(x[2]), 0)
    })

    insertions <- reactive({
        algLines()[seq(2, length(algLines()), 3)] |> gregexpr(pattern = '-+')
    })

    insertionsCollapse <- reactive({
        insertions() |> lapply(
            function(x) {
                x <- x - c(0, head(cumsum(attr(x, which = "match.length")), -1))
            }
        )
    })

    queryMd <- reactive({
        queryMd <- rep("", length(algLines()) / 3)
        for (i in seq_len(length(queryMd))) {
            query <- algLines()[3 * i]
            insertion <- insertions()[[i]]
            if (insertion[1] == -1) {
                queryMd[i] <- sprintf("<details><summary style=\"list-style-position: outside;\">%s</summary></details>", query)
                next
            }
            insertionLines <- arrangeInsertion(query, insertion, insertionsCollapse()[[i]])
            queryMd[i] <- "<details style=\"margin: 0, 0, 0, 0; padding: 0, 0, 0, 0;\"><summary style=\"list-style-position: outside;\">"
            if (insertion[1] == 1) {
                queryMd[i] <- sprintf("%s<span style=\"letter-spacing: -0.3em;\"> </span>", queryMd[i])
            } else {
                queryMd[i] <- sprintf("%s ", queryMd[i])
            }
            start <- 1
            for (j in seq_len(length(insertion))) {
                queryMd[i] <- sprintf("%s%s<span style=\"letter-spacing: -0.3em;\">%s|</span>", queryMd[i], substr(query, start, insertion[j] - 2), substr(query, insertion[j] - 1, insertion[j] - 1))
                start <- insertion[j] + attr(insertion, which = "match.length")[j]
            }
            queryMd[i] <- sprintf("%s%s</summary>%s</details>",
                queryMd[i],
                substr(query, start, nchar(query)),
                paste(insertionLines, collapse = '<br>')
            )
        }
        queryMd
    })

    refMd <- reactive({
        req(input$cut1)
        req(input$cut2)
        paste0(
            sprintf("<span style=\"background-color: lightgrey;\">%s", substr(reference(), 1, input$cut1 - 1)),
            sprintf('<span style=\"letter-spacing: -0.3em;\">%s|</span>', substr(reference(), input$cut1, input$cut1)),
            sprintf('<span style=\"color: red;\">%s</span>', substr(reference(), input$cut1 + 1, ref1Len())),
            sprintf('<span style=\"color: blue;\">%s</span>', substr(reference(), ref1Len() + 1, ref1Len() + input$cut2 - 1)),
            sprintf('<span style=\"color: blue; letter-spacing: -0.3em;\">%s</span><span style=\"letter-spacing: -0.3em;\">|</span>', substr(reference(), ref1Len() + input$cut2, ref1Len() + input$cut2)),
            sprintf("%s</span>", substr(reference(), ref1Len() + input$cut2 + 1, nchar(reference())))
        )
    })

    output$reference <- renderText({
        refMd()
    })

    output$queries <- renderText({
        paste(c(queryMd(), rep("<div>~</div>", 20)), collapse = "\n")
        
    })
}

# Run the app ----
shinyApp(ui = ui, server = server)
getIndelTypes <- function(algTibble) {
    algTibble |> mutate(
        indelType = factor(
            ifelse(
                insert & !delete,
                "insertion",
                ifelse(
                    delete & !insert,
                    "deletion",
                    ifelse(insert & delete,
                        "indel",
                        "WT"
                    )
                )
            ),
            levels = c("WT", "deletion", "insertion", "indel")
        ),
        count = count
    ) |> summarise(count = sum(count), .by = "indelType")
}

getIndelTypesEx <- function(algTibble) {
    algTibble |> mutate(
        indelType = factor(
            ifelse(templatedInsert & delete & nchar(randInsert),
                "full",
                ifelse(templatedInsert & delete & !nchar(randInsert),
                    "tempdel",
                    ifelse(templatedInsert & !delete & nchar(randInsert),
                        "temprand",
                        ifelse(templatedInsert & !delete & !nchar(randInsert),
                            "templated",
                            ifelse(!templatedInsert & delete & nchar(randInsert),
                                "randdel",
                                ifelse(!templatedInsert & delete & !nchar(randInsert),
                                    "deletion",
                                    ifelse(!templatedInsert & !delete & nchar(randInsert),
                                        "random",
                                        "WT"
                                    )
                                )
                            )
                        )
                    )
                )
            ),
            levels = c("WT", "deletion", "templated", "random", "temprand", "tempdel", "randdel", "full")
        ),
        count = count
    ) |> summarise(count = sum(count), .by = "indelType")
}

indelTypePiePlot <- function(indelTypeTibble, classifyTempFile) {
    fillColors <- RColorBrewer::brewer.pal(8, "Set1")
    if (nrow(indelTypeTibble) == 4) {
        fillColors <- fillColors[c(1, 2, 3, 6)]
    }
    ggFig <- indelTypeTibble |>
        mutate(percent = count / sum(count), perlabel = scales::percent(percent, accuracy = 0.01)) |>
        mutate(typeCount = sprintf("%s: %d", indelType, count)) |>
        mutate(
            typeCount = factor(
                typeCount,
                levels = sapply(
                    indelTypeTibble$indelType |> levels(),
                    function(level) typeCount[startsWith(typeCount, level)],
                    USE.NAMES = FALSE
                ) |> unlist()
            )
        ) |>
        ggplot(aes(1, percent, fill = typeCount, weight = count)) +
        geom_col() +
        geom_text(aes(label = perlabel), position = position_stack(vjust = 0.5), size = 5) +
        scale_x_discrete(name = NULL, breaks = NULL) +
        scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75), labels = scales::percent) +
        scale_fill_manual(values = fillColors) +
        coord_polar(theta = "y") +
        theme(text = element_text(size = 30))
    ggsave(paste0(classifyTempFile, ".pdf"), plot = ggFig)
    tags$iframe(src = paste0(sub("^www/", "", classifyTempFile), ".pdf"), height = "1200px", width = "100%")
}

indelTypeWafflePlot <- function(indelTypeTibble, classifyTempFile) {
    fillColors <- RColorBrewer::brewer.pal(8, "Set1")
    if (nrow(indelTypeTibble) == 4) {
        fillColors <- fillColors[c(1, 2, 3, 6)]
    }
    ggFig <- ggplot(indelTypeTibble, aes(fill = indelType, values = count)) +
        geom_waffle(n_rows = 10, color = "white", make_proportional = TRUE) +
        scale_fill_manual(limits = indelTypeTibble$indelType |> levels(), values = fillColors) +
        coord_equal() +
        theme(text = element_text(size = 30))
    ggsave(paste0(classifyTempFile, ".pdf"), plot = ggFig)
    tags$iframe(src = paste0(sub("^www/", "", classifyTempFile), ".pdf"), height = "1200px", width = "100%")
}
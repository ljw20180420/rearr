getIndelTypes <- function(inserts, deletes, counts) {
    tibble(
        indelType = factor(
            ifelse(
                inserts & !deletes,
                "insertion",
                ifelse(
                    deletes & !inserts,
                    "deletion",
                    ifelse(inserts & deletes,
                        "indel",
                        "WT"
                    )
                )
            ),
            levels = c("WT", "deletion", "insertion", "indel")
        ),
        count = counts
    ) |> summarise(count = sum(count), .by = "indelType")
}

getIndelTypesEx <- function(templatedInserts, deletes, randomInserts, counts) {
    tibble(
        indelType = factor(
            ifelse(templatedInserts & deletes & randomInserts,
                "full",
                ifelse(templatedInserts & deletes & !randomInserts,
                    "tempdel",
                    ifelse(templatedInserts & !deletes & randomInserts,
                        "temprand",
                        ifelse(templatedInserts & !deletes & !randomInserts,
                            "templated",
                            ifelse(!templatedInserts & deletes & randomInserts,
                                "randdel",
                                ifelse(!templatedInserts & deletes & !randomInserts,
                                    "deletion",
                                    ifelse(!templatedInserts & !deletes & randomInserts,
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
        count = counts
    ) |> summarise(count = sum(count), .by = "indelType")
}

indelTypePiePlot <- function(indelTypeTibble) {
    fillColors <- RColorBrewer::brewer.pal(8, "Set1")
    if (nrow(indelTypeTibble) == 4) {
        fillColors <- fillColors[c(1, 2, 3, 6)]
    }
    indelTypeTibble |>
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
}

indelTypeWafflePlot <- function(indelTypeTibble) {
    fillColors <- RColorBrewer::brewer.pal(8, "Set1")
    if (nrow(indelTypeTibble) == 4) {
        fillColors <- fillColors[c(1, 2, 3, 6)]
    }
    ggplot(indelTypeTibble, aes(fill = indelType, values = count)) +
    geom_waffle(n_rows = 10, color = "white", make_proportional = TRUE) +
    scale_fill_manual(limits = indelTypeTibble$indelType |> levels(), values = fillColors) +
    coord_equal() +
    theme(text = element_text(size = 30))
}
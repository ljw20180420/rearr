#' Classify editing output
#'
#' Classify the CRISPR editing output based on whether there is deletion and insertion.
#'
#' @param algTibble a tibble containing editing information
#' @return a summary tibble with the column indelType ("WT", "deletion", "insertion", "indel") and count (the count of indelType)
#' @export
#'
#' @examples
#' indelTypeTibble <- getIndelTypes(algTibble)
getIndelTypes <- function(algTibble) {
  algTibble |>
    dplyr::mutate(
      indelType = factor(
        ifelse(
          insert & !delete,
          "insertion",
          ifelse(
            delete & !insert,
            "deletion",
            ifelse(insert & delete, "indel", "WT")
          )
        ),
        levels = c("WT", "deletion", "insertion", "indel")
      ),
      count = count
    ) |>
    dplyr::summarise(count = sum(count), .by = "indelType")
}

#' Classify editing output with discrimination of templated and random insertions
#'
#' Classify the CRISPR editing output based on whether there is deletion, templated insertion and random insertion.
#'
#' @param algTibble a tibble containing editing information
#' @return a summary tibble with the column indelType ("WT", "deletion", "templated", "random", "temprand", "tempdel", "randdel", "full") and count (the count of indelType)
#' @export
#'
#' @examples
#' indelTypeTibbleEx <- getIndelTypesEx(algTibble)
getIndelTypesEx <- function(algTibble) {
  algTibble |>
    dplyr::mutate(
      indelType = factor(
        ifelse(
          templatedInsert & delete & nchar(randInsert),
          "full",
          ifelse(
            templatedInsert & delete & !nchar(randInsert),
            "tempdel",
            ifelse(
              templatedInsert & !delete & nchar(randInsert),
              "temprand",
              ifelse(
                templatedInsert & !delete & !nchar(randInsert),
                "templated",
                ifelse(
                  !templatedInsert & delete & nchar(randInsert),
                  "randdel",
                  ifelse(
                    !templatedInsert & delete & !nchar(randInsert),
                    "deletion",
                    ifelse(
                      !templatedInsert & !delete & nchar(randInsert),
                      "random",
                      "WT"
                    )
                  )
                )
              )
            )
          )
        ),
        levels = c(
          "WT",
          "deletion",
          "templated",
          "random",
          "temprand",
          "tempdel",
          "randdel",
          "full"
        )
      ),
      count = count
    ) |>
    dplyr::summarise(count = sum(count), .by = "indelType")
}

#' Pie plot of editing classfication
#'
#' Visualize the classification of editing output by pie plot.
#'
#' @param indelTypeTibble the tibble with classification and counts
#' @param classifyTempFile the pdf tempfile
#' @return the html iframe to be displayed by shiny:renderUI
#' @export
#'
#' @examples
#' iframeTag <- indelTypePiePlot(indelTypeTibble, classifyTempFile)
indelTypePiePlot <- function(indelTypeTibble, classifyTempFile) {
  fillColors <- RColorBrewer::brewer.pal(8, "Set1")
  if (nrow(indelTypeTibble) == 4) {
    fillColors <- fillColors[c(1, 2, 3, 6)]
  }
  ggFig <- indelTypeTibble |>
    dplyr::mutate(
      percent = count / sum(count),
      perlabel = scales::percent(percent, accuracy = 0.01)
    ) |>
    dplyr::mutate(typeCount = sprintf("%s: %d", indelType, count)) |>
    dplyr::mutate(
      typeCount = factor(
        typeCount,
        levels = sapply(
          indelTypeTibble$indelType |> levels(),
          function(level) typeCount[startsWith(typeCount, level)],
          USE.NAMES = FALSE
        ) |>
          unlist()
      )
    ) |>
    ggplot2::ggplot(ggplot2::aes(
      1,
      percent,
      fill = typeCount,
      weight = count
    )) +
    ggplot2::geom_col() +
    ggplot2::geom_text(
      ggplot2::aes(label = perlabel),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 5
    ) +
    ggplot2::scale_x_discrete(name = NULL, breaks = NULL) +
    ggplot2::scale_y_continuous(
      breaks = c(0, 0.25, 0.5, 0.75),
      labels = scales::percent
    ) +
    ggplot2::scale_fill_manual(values = fillColors) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::theme(text = ggplot2::element_text(size = 30))
  ggplot2::ggsave(
    classifyTempFile,
    plot = ggFig,
    height = 3600,
    width = 3600,
    unit = "px"
  )
  htmltools::tags$iframe(
    src = sub("^www/", "", classifyTempFile),
    height = "1200px",
    width = "100%"
  )
}

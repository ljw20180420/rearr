discreteDistribution <- function(algTibble, distriTempFile) {
  ggFig <- algTibble |>
    tidyr::pivot_longer(
      cols = -"count",
      names_to = "type",
      values_to = "value"
    ) |>
    dplyr::mutate(type = factor(type, levels = colnames(algTibble))) |>
    ggplot2::ggplot(ggplot2::aes(x = value, weight = count, fill = type)) +
    ggplot2::geom_bar(position = ggplot2::position_dodge()) +
    ggplot2::scale_y_continuous(name = "count", expand = c(0, 0))
  ggplot2::ggsave(distriTempFile, plot = ggFig)
  htmltools::tags$iframe(
    src = sub("^www/", "", distriTempFile),
    height = "1200px",
    width = "100%"
  )
}

continuousDistribution <- function(algTibble, distriTempFile) {
  ggFig <- algTibble |>
    tidyr::pivot_longer(
      cols = -"count",
      names_to = "type",
      values_to = "value"
    ) |>
    dplyr::mutate(type = factor(type, levels = colnames(algTibble))) |>
    ggplot2::ggplot(ggplot2::aes(value, color = type, weight = count)) +
    ggplot2::geom_density(ggplot2::aes(y = ggplot2::after_stat(density))) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(name = "density", expand = c(0, 0))
  ggplot2::ggsave(distriTempFile, plot = ggFig)
  htmltools::tags$iframe(
    src = sub("^www/", "", distriTempFile),
    height = "1200px",
    width = "100%"
  )
}

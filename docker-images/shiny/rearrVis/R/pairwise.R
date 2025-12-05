pairwisePlot <- function(pairTibble, xscale, yscale, method, span) {
  ggFig <- pairTibble |>
    ggplot2::ggplot(ggplot2::aes(
      x = .data[[colnames(pairTibble)[1]]],
      y = .data[[colnames(pairTibble)[2]]]
    )) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = method, span = span) +
    ggplot2::scale_x_continuous(expand = c(0, 0), transform = xscale) +
    ggplot2::scale_y_continuous(expand = c(0, 0), transform = yscale)
}

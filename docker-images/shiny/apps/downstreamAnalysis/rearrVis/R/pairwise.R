#' Pairwise plot for two columns
#'
#' Scatter plot with smooth fit for two columns.
#'
#' @param pairTibble tibble with two columns to apply pairwise plot
#' @param xscale transform of x axis
#' @param yscale transform of y axis
#' @param method method for smooth fit
#' @param span span for smooth fit
#' @param pairwiseTempFile the pdf temp file
#' @return the html iframe tag to be displayed in shiny::renderUI
#' @export
#'
#' @examples
#' # example code
pairwisePlot <- function(
  pairTibble,
  xscale,
  yscale,
  method,
  span,
  pairwiseTempFile
) {
  ggFig <- pairTibble |>
    ggplot2::ggplot(ggplot2::aes(
      x = .data[[colnames(pairTibble)[1]]],
      y = .data[[colnames(pairTibble)[2]]]
    )) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = method, formula = "y ~ x", span = span) +
    ggplot2::scale_x_continuous(expand = c(0, 0), transform = xscale) +
    ggplot2::scale_y_continuous(expand = c(0, 0), transform = yscale)
  suppressMessages(
    ggplot2::ggsave(
      pairwiseTempFile,
      ggFig,
      height = 3600,
      width = 3600,
      units = "px"
    )
  )
  htmltools::tags$iframe(
    src = sub("^www/", "", pairwiseTempFile),
    height = "1200px",
    width = "100%"
  )
}

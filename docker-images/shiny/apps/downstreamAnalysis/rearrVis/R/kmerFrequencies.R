#' Kmer enrichment of target
#'
#' Stack the target and nontarget editing output count for each kmer from sgRNA.
#'
#' @param kmerTibble tibble with kmer count for target and nontarget editing output
#' @param kmerPdfTempFile the pdf temp file
#' @return the html iframe tag to be displayed in shiny::renderUI
#' @export
#'
#' @examples
#' iframeTag <- plotKmerFrequencies(kmerTibble, kmerPdfTempFile)
plotKmerFrequencies <- function(kmerTibble, kmerPdfTempFile) {
  ggFig <- kmerTibble |>
    ggplot2::ggplot(ggplot2::aes(x = kmer, y = count, fill = target)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 45)) +
    ggplot2::scale_y_continuous(expand = c(0, 0))
  ggplot2::ggsave(
    kmerPdfTempFile,
    plot = ggFig,
    height = 3600,
    width = 3600,
    units = "px"
  )
  htmltools::tags$iframe(
    src = sub("^www/", "", kmerPdfTempFile),
    height = "1200px",
    width = "100%"
  )
}

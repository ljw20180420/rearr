plotKmerFrequencies <- function(kmerTibble, kmerPdfTempFile) {
  ggFig <- kmerTibble |>
    ggplot2::ggplot(ggplot2::aes(kmer, count, fill = target)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 45)) +
    ggplot2::scale_y_continuous(expand = c(0, 0))
  ggplot2::ggsave(kmerPdfTempFile, plot = ggFig)
  htmltools::tags$iframe(
    src = sub("^www/", "", kmerPdfTempFile),
    height = "1200px",
    width = "100%"
  )
}

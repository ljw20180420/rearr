plotKmerFrequencies <- function(kmerTibble, kmerPdfTempFile) {
    ggFig <- kmerTibble |>
        ggplot(aes(kmer, count, fill = target)) +
        geom_col(position = "stack") +
        scale_x_discrete(guide = guide_axis(angle = 45)) +
        scale_y_continuous(expand = c(0, 0))
    ggsave(kmerPdfTempFile, plot = ggFig)
    tags$iframe(src=sub("^www/", "", kmerPdfTempFile), height="1200px", width = "100%")
}
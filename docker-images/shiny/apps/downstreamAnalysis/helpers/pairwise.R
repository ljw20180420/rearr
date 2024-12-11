pairwisePlot <- function(pairTibble, xscale, yscale, method, span) {
    ggFig <- pairTibble |>
        ggplot(aes(x = .data[[colnames(pairTibble)[1]]], y = .data[[colnames(pairTibble)[2]]])) +
        geom_point() +
        geom_smooth(method = method, span = span) +
        scale_x_continuous(expand = c(0, 0), transform = xscale) +
        scale_y_continuous(expand = c(0, 0), transform = yscale)
}
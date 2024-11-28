discreteDistribution <- function(algTibble, distriTempFile) {
    ggFig <- algTibble |>
        pivot_longer(cols = -"count", names_to = "type", values_to = "value") |>
        mutate(type = factor(type, levels = colnames(algTibble))) |>
        ggplot(aes(x = value, weight = count, fill = type)) +
        geom_bar(position = position_dodge()) +
        scale_y_continuous(name = "count", expand = c(0, 0))
    ggsave(paste0(distriTempFile, ".pdf"), plot = ggFig)
    tags$iframe(src = paste0(sub("^www/", "", distriTempFile), ".pdf"), height = "1200px", width = "100%")
}

continuousDistribution <- function(algTibble, distriTempFile) {
    ggFig <- algTibble |>
        pivot_longer(cols = -"count", names_to = "type", values_to = "value") |>
        mutate(type = factor(type, levels = colnames(algTibble))) |>
        ggplot(aes(value, color = type, weight = count)) +
        geom_density(aes(y = after_stat(density))) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(name = "density", expand = c(0, 0))
    ggsave(paste0(distriTempFile, ".pdf"), plot = ggFig)
    tags$iframe(src = paste0(sub("^www/", "", distriTempFile), ".pdf"), height = "1200px", width = "100%")
}
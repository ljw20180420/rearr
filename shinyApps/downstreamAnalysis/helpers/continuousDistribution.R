discreteDistribution <- function(distriTibble) {
    distriTibble |>
        pivot_longer(cols = -"count", names_to = "type", values_to = "value") |>
        mutate(type = factor(type, levels = colnames(distriTibble))) |>
        ggplot(aes(x = value, weight = count, fill = type)) +
        geom_bar(position = position_dodge()) +
        scale_y_continuous(name = "count", expand = c(0, 0))
}

continuousDistribution <- function(distriTibble) {
    distriTibble |>
        pivot_longer(cols = -"count", names_to = "type", values_to = "value") |>
        mutate(type = factor(type, levels = colnames(distriTibble))) |>
        ggplot(aes(value, color = type, weight = count)) +
        geom_density(aes(y = after_stat(density))) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(name = "density", expand = c(0, 0))
}

getDistriTibble <- function(counts, scores, cuts1, cuts2, ref1Lens, ref2Lens, ref1ends, ref2starts, upInserts, downInserts, randomInserts, templatedInserts, inserts, upDeletes, downDeletes, deletes) {
    distriTibble <- tibble(
        count = counts,
        "alignment score" = scores,
        cut1 = cuts1,
        cut2 = cuts2,
        "ref1 length" = ref1Lens,
        "ref2 length" = ref2Lens,
        "alignment end in ref1" = ref1ends - cuts1,
        "alignment start in ref2" = ref2starts - cuts2,
        "upstream insertion" = upInserts,
        "downstream insertion" = downInserts,
        "random insertion" = randomInserts,
        "templated insertion" = templatedInserts,
        "insertion" = inserts,
        "upstream deletion" = upDeletes,
        "downstream deletion" = downDeletes,
        "deletion" = deletes
    )
    return(distriTibble)
}
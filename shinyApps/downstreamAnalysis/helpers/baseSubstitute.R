countBaseSubstitute <- function(algLines, counts) {
    refBases <- algLines[seq(2, length(algLines), 3)] |> toupper() |> strsplit("") |> unlist()
    queryBases <- algLines[seq(3, length(algLines), 3)] |> strsplit("") |> unlist()
    refLens <- vapply(algLines[seq(2, length(algLines), 3)], nchar, 0, USE.NAMES = FALSE)
    levels <- paste(rep(c("-", "A", "C", "G", "T"), times = 5), rep(c("-", "A", "C", "G", "T"), each = 5), sep = ">")
    return(
        tibble(
            sub = paste(refBases, queryBases, sep = ">"),
            count = rep(counts, times = refLens)
        ) |>
        summarise(count = sum(count), .by = "sub") |>
        mutate(sub = factor(sub, levels = levels))
    )
}
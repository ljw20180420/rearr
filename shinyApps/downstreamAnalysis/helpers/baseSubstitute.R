countBaseSubstitute <- function(algTibble) {
    levels <- paste(rep(c("-", "A", "C", "G", "T"), times = 5), rep(c("-", "A", "C", "G", "T"), each = 5), sep = ">")
    tibble(
        sub = lapply(seq_len(nrow(algTibble)), function(i) {
            paste(
                algTibble$refLine[i] |> toupper() |> str_split_1(""),
                algTibble$queryLine[i] |> str_split_1(""),
                sep = ">"
            )
        }),
        count = algTibble$count
    ) |> unnest(sub) |> summarise(count = sum(count), .by = "sub") |> mutate(sub = factor(sub, levels = levels))
}
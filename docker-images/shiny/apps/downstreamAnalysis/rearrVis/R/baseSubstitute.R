#' Count base substitution of query in reference.
#'
#' Mutation is like "A>C". Insertion is like "->A". Deletion is like "A>-".
#'
#' @param algTibble tibble containing gaggped reference and query in the two-line alignment.
#' @return tibble with sub labels like X>X and corresponding counts.
#' @export
#'
#' @examples
#' tb <- countBaseSubstitute(algTibble)
countBaseSubstitute <- function(algTibble) {
  levels <- paste(
    rep(c("-", "A", "C", "G", "T"), times = 5),
    rep(c("-", "A", "C", "G", "T"), each = 5),
    sep = ">"
  )
  tibble::tibble(
    sub = lapply(seq_len(nrow(algTibble)), function(i) {
      paste(
        algTibble$refLine[i] |> toupper() |> stringr::str_split_1(""),
        algTibble$queryLine[i] |> stringr::str_split_1(""),
        sep = ">"
      )
    }),
    count = algTibble$count
  ) |>
    tidyr::unnest(sub) |>
    dplyr::summarise(count = sum(count), .by = "sub") |>
    dplyr::mutate(sub = factor(sub, levels = levels))
}

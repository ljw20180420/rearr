extendToSameLength <- function(seqs) {
    maxLen <- max(nchar(seqs))
    return(
        seqs |> str_pad(width = maxLen, side = 'right', pad = "-") |>
        strsplit("") |> unlist() |> matrix(ncol = maxLen, byrow = TRUE)
    )
}

extendToAlignCut <- function(seqs, cuts) {
    leftMax <- max(cuts)
    rightMax <- max(nchar(seqs) - cuts)
    return(
        paste0(
            substr(seqs, 1, cuts) |> str_pad(width = leftMax, side = 'left', pad = "-"),
            substr(seqs, cuts + 1, nchar(seqs)) |> str_pad(width = rightMax, side = 'right', pad = "-")
        ) |> strsplit("") |> unlist() |> matrix(ncol = leftMax + rightMax, byrow = TRUE)
    )
}

vectorToStringVector <- function(vec, strLens) {
    split(vec, rep.int(seq_along(strLens), strLens)) |> vapply(function(x) paste0(x, collapse = ""), "", USE.NAMES = FALSE)
}

posMatrixToTibble <- function(mat, cut) {
    tibble(
        count = c(mat),
        pos = rep(seq(ncol(mat)) - 0.5 - cut, each = nrow(mat)),
        type = factor(rep(rownames(mat), times = ncol(mat)), levels = rownames(mat))
    )
}

drawPositionalStatic <- function(inputTibble, insertCount) {
    inputTibble |>
        ggplot(aes(pos, count)) +
        geom_col(aes(fill = type)) +
        geom_step(aes(pos, count, color = "black"), data = insertCount[2:(nrow(insertCount) - 1),], direction = "mid") +
        scale_x_continuous(name = "position relative to cut1", expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_color_identity(name = NULL, guide = guide_legend(), labels = "insertion")
}
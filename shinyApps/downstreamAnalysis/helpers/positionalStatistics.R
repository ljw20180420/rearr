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
        geom_col(aes(fill = type), width = 1) +
        geom_step(aes(pos, count, color = "black"), data = insertCount[2:(nrow(insertCount) - 1),], direction = "mid") +
        scale_x_continuous(name = "position relative to cut1", expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_color_identity(name = NULL, guide = guide_legend(), labels = "insertion")
}

getPositionalBaseFreq <- function(queryMat, counts) {
    baseFreq <- rbind(
        colSums((queryMat == "-") * counts),
        colSums((queryMat == "A") * counts),
        colSums((queryMat == "C") * counts),
        colSums((queryMat == "G") * counts),
        colSums((queryMat == "T") * counts)
    )
    rownames(baseFreq) <- c("-", "A", "C", "G", "T")
    return(baseFreq)
}

getPositionalMSDTibble <- function(queryMat, refMat, counts, cut) {
    refMat <- toupper(refMat)
    MSDFreq <- rbind(
        colSums((queryMat == '-') * counts),
        colSums((queryMat != refMat & queryMat != '-') * counts),
        colSums((queryMat == refMat) * counts)
    )
    rownames(MSDFreq) <- c("delete", "SNP", "match")
    MSDFreq |> posMatrixToTibble(cut)
}

calInsertionCount <- function(refList, cuts, refLens) {
    insertList <- vector("list", length(refList))
    for (i in seq_len(length(refList))) {
        mask <- refList[[i]] != "-"
        insertList[[i]] <- cumsum(mask)[!mask] - cuts[i]
    }
    histCount <- insertList |> unlist() |> table()
    fullCount <- rep(0, max(cuts) + max(refLens - cuts) + 1)
    fullCount[as.integer(names(histCount)) + max(cuts) + 1] <- histCount
    tibble(
        count = fullCount,
        pos = seq(-max(cuts), max(refLens - cuts))
    )
}

drawPositionalReads <- function(queryMat, cut) {
    queryMat |> melt() |>
        ggplot(aes(x = Var2 - cut, y = Var1)) + 
        geom_raster(aes(fill=value), hjust = 0, vjust = 0) + 
        scale_fill_manual(values = c("-" = "white", "A" = "darkgreen", "C" = "blue", "G" = "gold", "T" = "red")) +
        scale_x_continuous(name = "position relative to cut1", expand = c(0, 0)) +
        scale_y_continuous(name = "reads", expand = c(0, 0))
}

drawPositionalSnps <- function(queryMat, refMat, cut) {
    refMat <- toupper(refMat)
    queryMat[queryMat != refMat & queryMat != "-"] = "S"
    queryMat[queryMat == "-"] = "D"
    queryMat[queryMat == refMat] = "M"
    queryMat |> melt() |>
        ggplot(aes(x = Var2 - cut, y = Var1)) + 
        geom_raster(aes(fill=value), hjust = 0, vjust = 0) + 
        scale_fill_manual(values = c("D" = "white", "M" = "darkgreen", "S" = "red")) +
        scale_x_continuous(name = "position relative to cut1", expand = c(0, 0)) +
        scale_y_continuous(name = "reads", expand = c(0, 0))
}

drawPositionalLogo <- function(baseFreq, method, namespace) {
    ggplot() +
    geom_logo(data = baseFreq, method = method, namespace = namespace) +
    scale_x_continuous(breaks = NULL) +
    theme_logo()
}
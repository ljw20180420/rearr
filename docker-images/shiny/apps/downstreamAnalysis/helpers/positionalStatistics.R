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

drawPositionalStatic <- function(inputTibble, insertCount, posBaseRefTempFile) {
    ggFig <- inputTibble |>
        ggplot(aes(pos, count)) +
        geom_col(aes(fill = type), width = 1) +
        geom_step(aes(pos, count, color = "black"), data = insertCount[2:(nrow(insertCount) - 1),], direction = "mid") +
        scale_x_continuous(name = "position relative to cut", expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_color_identity(name = NULL, guide = guide_legend(), labels = "insertion")
    ggsave(posBaseRefTempFile, plot = ggFig, height = 1200, width = 3600, unit = "px")
    tags$iframe(src = sub("^www/", "", posBaseRefTempFile), height = "600px", width = "100%")
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

calInsertionCount <- function(refList, cuts, maxCutDown) {
    maxCut <- max(cuts)
    insertList <- vector("list", length(refList))
    for (i in seq_len(length(refList))) {
        mask <- refList[[i]] != "-"
        insertList[[i]] <- cumsum(mask)[!mask] - cuts[i]
    }
    histCount <- insertList |> unlist() |> table()
    fullCount <- rep(0, maxCut + maxCutDown + 1)
    fullCount[as.integer(names(histCount)) + maxCut + 1] <- histCount
    tibble(
        count = fullCount,
        pos = seq(-maxCut, maxCutDown)
    )
}

resolutionMatrix <- function(mat, counts, resolution) {
    totalCount <- sum(counts)
    resolutionCount <- round(seq(1, totalCount + 1, length.out = resolution + 1))
    resolutionMatrix <- matrix(0, resolution, ncol(mat))
    cumCount <- 0
    j <- 1
    for (i in seq_len(length(counts))) {
        cumCount <- cumCount + counts[i]
        if (cumCount < resolutionCount[j + 1]) {
            resolutionMatrix[j,] <- resolutionMatrix[j,] + mat[i,] * counts[i]
        } else {
            resolutionMatrix[j,] <- resolutionMatrix[j,] + mat[i,] * (resolutionCount[j + 1] - 1 + counts[i] - cumCount)
            j <- j + 1
            while (resolutionCount[j + 1] <= cumCount) {
                resolutionMatrix[j,] <- mat[i,] * (resolutionCount[j + 1] - resolutionCount[j])
                j <- j + 1
            }
            resolutionMatrix[j,] <- mat[i,] * (cumCount - resolutionCount[j] + 1)
        }
    }
    return(resolutionMatrix)
}

getPositionalReads <- function(mat, counts, resolution, cut) {
    resolutionMat <- resolutionMatrix(mat, counts, resolution)
    resolutionCount <- round(seq(1, sum(counts) + 1, length.out = resolution + 1))
    vertCent <- (resolutionCount[-length(resolutionCount)] + resolutionCount[-1]) / 2 - 1
    resolutionMat |> melt() |> mutate(x = Var2 - cut - 0.5, y = vertCent[Var1], value = value, height = resolutionCount[Var1 + 1] - resolutionCount[Var1]) |> select(x, y, value, height)
}

drawPositionalReads <- function(tibb, maxCut, maxCutdown, posBaseRefTempFile) {
    ggFig <- ggplot(tibb, aes(x = x, y = y, fill = value, height = height)) + 
        geom_tile(width = 1) + 
        scale_fill_gradient(low = "white", high = "red") +
        scale_x_continuous(name = "position relative to cut", limits = c(-maxCut, maxCutdown), expand = c(0, 0)) +
        scale_y_continuous(name = "reads", expand = c(0, 0))
    ggsave(posBaseRefTempFile, plot = ggFig, height = 1200, width = 3600, unit = "px")
    tags$iframe(src = sub("^www/", "", posBaseRefTempFile), height = "600px", width = "100%")
}

drawPositionalLogo <- function(baseFreq, method, namespace, posBaseRefTempFile) {
    ggFig <- ggplot() +
        geom_logo(data = baseFreq, method = method, namespace = namespace) +
        scale_x_continuous(breaks = NULL) +
        theme_logo()
    ggsave(posBaseRefTempFile, plot = ggFig, height = 1200, width = 3600, unit = "px")
    tags$iframe(src = sub("^www/", "", posBaseRefTempFile), height = "600px", width = "100%")
}
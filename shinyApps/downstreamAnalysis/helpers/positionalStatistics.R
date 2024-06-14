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
    ggsave(paste0(posBaseRefTempFile, ".pdf"), plot = ggFig, height = 1200, width = 3600, unit = "px")
    tags$iframe(src = paste0(sub("^www/", "", posBaseRefTempFile), ".pdf"), height = "600px", width = "100%")
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

downSampleMatrix <- function(mat, counts, thresCount) {
    cls <- 1
    cumCount <- 0
    clses <- rep(0, length(counts))
    for (i in seq_len(length(counts))) {
        clses[i] <- cls
        cumCount <- cumCount + counts[i]
        if (cumCount > thresCount) {
            cls <- cls + 1
            cumCount <- 0
        }
    }
    mat <- as_tibble(mat)
    mat |> mutate(count = counts, cls = clses) |> summarise(across(colnames(mat), ~ .x[which.max(count)]), count = sum(count), .by = "cls") |> select(-"cls")
}

getPositionalReads <- function(queryMat, counts, thresCount, cut) {
    queryTibbleDown <- downSampleMatrix(queryMat, counts, thresCount)
    cumCounts <- cumsum(queryTibbleDown$count)
    vertCent <- c(cumCounts[1] / 2, (cumCounts[-length(cumCounts)] + cumCounts[-1]) / 2)
    queryTibbleDown |> select(-"count") |> as.matrix() |> `colnames<-`(seq_len(ncol(queryTibbleDown) - 1)) |> melt() |> filter(value != "-") |> mutate(x = Var2 - cut - 0.5, y = vertCent[Var1], value = value, height = queryTibbleDown$count[Var1]) |> select(x, y, value, height)
}

drawPositionalReads <- function(readTibble, maxCut, maxCutdown, posBaseRefTempFile) {
    ggFig <- ggplot(readTibble, aes(x = x, y = y, fill = value, height = height)) + 
        geom_tile(width = 1) + 
        scale_fill_manual(values = c("A" = "darkgreen", "C" = "blue", "G" = "gold", "T" = "red")) +
        scale_x_continuous(name = "position relative to cut", limits = c(-maxCut, maxCutdown), expand = c(0, 0)) +
        scale_y_continuous(name = "reads", expand = c(0, 0))
    ggsave(paste0(posBaseRefTempFile, ".pdf"), plot = ggFig, height = 1200, width = 3600, unit = "px")
    tags$iframe(src = paste0(sub("^www/", "", posBaseRefTempFile), ".pdf"), height = "600px", width = "100%")
}

drawPositionalSnps <- function(snpTibble, maxCut, maxCutdown, posBaseRefTempFile) {
    ggFig <- ggplot(snpTibble, aes(x = x, y = y, fill = value, height = height)) + 
        geom_tile(width = 1) + 
        scale_fill_manual(values = c("M" = "darkgreen", "S" = "red")) +
        scale_x_continuous(name = "position relative to cut", limits = c(-maxCut, maxCutdown), expand = c(0, 0)) +
        scale_y_continuous(name = "reads", expand = c(0, 0))
    ggsave(paste0(posBaseRefTempFile, ".pdf"), plot = ggFig, height = 1200, width = 3600, unit = "px")
    tags$iframe(src = paste0(sub("^www/", "", posBaseRefTempFile), ".pdf"), height = "600px", width = "100%")
}

drawPositionalLogo <- function(baseFreq, method, namespace, posBaseRefTempFile) {
    ggFig <- ggplot() +
        geom_logo(data = baseFreq, method = method, namespace = namespace) +
        scale_x_continuous(breaks = NULL) +
        theme_logo()
    ggsave(paste0(posBaseRefTempFile, ".pdf"), plot = ggFig, height = 1200, width = 3600, unit = "px")
    tags$iframe(src = paste0(sub("^www/", "", posBaseRefTempFile), ".pdf"), height = "600px", width = "100%")
}
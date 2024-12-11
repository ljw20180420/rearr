getMicroHomologyTibble <- function(ref1, ref2, cut1, cut2) {
    ref1vec <- str_split(toupper(ref1), "")[[1]]
    ref2vec <- str_split(toupper(ref2), "")[[1]]
    mh_matrix <- matrix(as.integer(rep(ref1vec, time = length(ref2vec)) == rep(ref2vec, each = length(ref1vec))), nrow = length(ref1vec))
    cls <- 1
    for (i in seq(2, nrow(mh_matrix))) {
        for (j in seq(2, ncol(mh_matrix))) {
            if (mh_matrix[i, j] > 0) {
                mh_matrix[i, j] <- mh_matrix[i - 1, j - 1] + mh_matrix[i, j]
            } else if (mh_matrix[i - 1, j - 1] > 0) {
                for (l in seq_len(mh_matrix[i - 1, j - 1])) {
                    mh_matrix[i - l, j - l] <- cls
                }
                cls <- cls + 1
            }
        }
    }
    for (i in seq(1, nrow(mh_matrix))) {
        if (mh_matrix[i, ncol(mh_matrix)] > 0) {
            for (l in seq_len(mh_matrix[i, ncol(mh_matrix)])) {
                mh_matrix[i - l + 1, ncol(mh_matrix) - l + 1] <- cls
            }
            cls <- cls + 1
        }
    }
    for (j in seq(1, ncol(mh_matrix) - 1)) {
        if (mh_matrix[nrow(mh_matrix), j] > 0) {
            for (l in seq_len(mh_matrix[nrow(mh_matrix), j])) {
                mh_matrix[nrow(mh_matrix) - l + 1, j - l + 1] <- cls
            }
            cls <- cls + 1
        }
    }
    return(
        mh_matrix |> melt() |> `colnames<-`(c("pos1", "pos2", "cls")) |> filter(cls > 0) |> summarise(pos1low = min(pos1) - 1 - cut1, pos1up = max(pos1) - cut1, shift = pos1up - max(pos2) + cut2, .by = "cls")
    )
}

drawMicroHomologyHeatmap <- function(mhTibbleSub, refEnd1Start2TibbleMicro, maxCut1, maxCut2, maxCut1down, maxCut2down, mode, mhMatrixTempFile) {
    if (mode == "separate") {
        refEnd1Start2TibbleMicro <- bind_rows(
            refEnd1Start2TibbleMicro |> filter(cls == 0),
            refEnd1Start2TibbleMicro |> filter(cls != 0) |> group_by(cls) |> mutate(count = count / n()) |> ungroup()
        )
    }
    mhPosTibble <- mhTibbleSub |> mutate(nacol = NA) |> pivot_longer(cols = c("pos1low", "pos1up", "nacol"), values_to = "mhPos1") |> mutate(mhPos2 = mhPos1 - shift) |> select(mhPos1, mhPos2)
    ggFig <- ggplot(refEnd1Start2TibbleMicro) +
        geom_tile(aes(x = pos2, y = pos1, fill = log10(count + 1)), height = 1, width = 1) +
        geom_path(aes(x = mhPos2, y = mhPos1), data = mhPosTibble) +
        scale_x_continuous(limits = c(-maxCut2 - 1, maxCut2down + 1), expand=c(0, 0)) +
        scale_y_continuous(limits = c(-maxCut1 - 1, maxCut1down + 1), expand=c(0, 0)) +
        scale_fill_gradientn(limits = c(0, NA), colors = c("blue", "white", "red")) +
        scale_size_area(max_size = 2) +
        coord_equal(ratio = 1)
    ggsave(mhMatrixTempFile, plot = ggFig)
    tags$iframe(src = sub("^www/", "", mhMatrixTempFile), height = "1200px", width = "100%")
}

getRefEnd1Start2Tibble <- function(algTibble, microRefId) {
    algTibble |> mutate(
        pos1 = ref1End - cut1,
        pos2 = ref2Start - cut2,
        refId = refId,
        count = count,
    ) |> filter(refId == microRefId) |> summarise(count = sum(count), .by = c("pos1", "pos2")) |> mutate(shift = pos1 - pos2)
}

getRefEnd1Start2TibbleMicro <- function(refEnd1Start2Tibble, mhTibbleSub) {
    joinTibble <- refEnd1Start2Tibble |> left_join(mhTibbleSub, by = "shift", relationship = "many-to-many")
    outRangeTibble <- joinTibble |> summarise(inRange = any(pos1 >= pos1low & pos1 <= pos1up), count = first(count), .by = c("pos1", "pos2")) |> filter(is.na(inRange) | !inRange) |> mutate(inRange = NULL, cls = 0)
    inRangeTibble <- joinTibble |> filter(pos1 >= pos1low, pos1 <= pos1up)
    if (nrow(inRangeTibble) == 0) {
        return(outRangeTibble)
    }
    inRangeTibble |> rowwise() |> mutate(pos1 = list(seq(pos1low, pos1up))) |> ungroup() |> unnest(pos1) |> mutate(pos2 = pos1 - shift) |> select(pos1, pos2, count, cls) |> bind_rows(outRangeTibble)
}
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
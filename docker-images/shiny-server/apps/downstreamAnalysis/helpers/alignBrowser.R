arrangeInsertion <- function(query, insertion, insertionCollapse) {
    insertionLines <- c("")
    for (i in seq_len(length(insertionCollapse))) {
        insertionAdded <- FALSE
        insertionSeg <- substr(query, insertion[i], insertion[i] + attr(insertion, which = "match.length")[i] - 1)
        for (j in seq_len(length(insertionLines))) {
            if (nchar(insertionLines[j]) < insertionCollapse[i]) {
                insertionLines[j] <- sprintf("%s%s%s ",
                    insertionLines[j],
                    strrep(" ", insertionCollapse[i] - nchar(insertionLines[j]) - 1),
                    insertionSeg
                )
                insertionAdded <- TRUE
                break
            }
        }
        if (!insertionAdded)
        {
            insertionLines <- c(
                insertionLines,
                sprintf("%s%s",
                    strrep(" ", insertionCollapse[i] - 1),
                    insertionSeg
                ) 
            )
        }
    }
    # add a = at the beginning to fit the start bar trick of query
    insertionLines <- paste0(sprintf("%10s", "#"), insertionLines)
    return(insertionLines)
}

snpHighlight <- function(refSeg, querySeg) {
    refArray <- str_split_1(toupper(refSeg), "")
    queryArray <- str_split_1(toupper(querySeg), "")
    diffs <- (refArray != queryArray & queryArray != "-") |> diff()
    diffs <- c(0, diffs, 0)
    boundaries <- diffs |> as.logical() |> which()
    if (length(boundaries) == 0) {
        return(querySeg)
    }
    mdStr <- ""
    start <- 1
    for (i in seq_len(length(boundaries) / 2)) {
        mdStr <- sprintf("%s%s<span style=\"color: red;\">%s</span>",
            mdStr,
            substr(querySeg, start, boundaries[2 * i - 1] - 1),
            substr(querySeg, boundaries[2 * i - 1], boundaries[2 * i] - 1)
        )
        start <- boundaries[2 * i]
    }
    mdStr <- paste0(mdStr, substr(querySeg, start, nchar(querySeg)))
}

getMarkdownFromAlign <- function(algTibble) {
    allMd <- rep("", nrow(algTibble) * 2)
    for (i in seq_len(nrow(algTibble))) {
        allMd[2 * i - 1] <- paste0(
            sprintf("<span style=\"background-color: lightgrey;\">%10s%s", algTibble$refId[i], substr(algTibble$refNoGap[i], 1, algTibble$cut1[i] - 1)),
            sprintf('<span style=\"letter-spacing: -0.3em;\">%s|</span>', substr(algTibble$refNoGap[i], algTibble$cut1[i], algTibble$cut1[i])),
            sprintf('<span style=\"color: red;\">%s</span>', substr(algTibble$refNoGap[i], algTibble$cut1[i] + 1, algTibble$ref1Len[i])),
            sprintf('<span style=\"color: blue;\">%s</span>', substr(algTibble$refNoGap[i], algTibble$ref1Len[i] + 1, algTibble$ref1Len[i] + algTibble$cut2[i] - 1)),
            sprintf('<span style=\"color: blue; letter-spacing: -0.3em;\">%s</span><span style=\"letter-spacing: -0.3em;\">|</span>', substr(algTibble$refNoGap[i], algTibble$ref1Len[i] + algTibble$cut2[i], algTibble$ref1Len[i] + algTibble$cut2[i])),
            sprintf("%s</span>", substr(algTibble$refNoGap[i], algTibble$ref1Len[i] + algTibble$cut2[i] + 1, nchar(algTibble$refNoGap[i])))
        )
        insertion <- algTibble$refLine[i] |> gregexpr(pattern = '-+') |> _[[1]]
        if (insertion[1] == -1) {
            queryMd <- sprintf("<details><summary style=\"list-style-position: outside;\">%10d%s</summary></details>", algTibble$count[i], algTibble$queryLine[i])
            allMd[2 * i] <- queryMd
            next
        }
        insertionCollapse <- insertion - c(0, head(cumsum(attr(insertion, which = "match.length")), -1))
        insertionLines <- arrangeInsertion(algTibble$queryLine[i], insertion, insertionCollapse)
        queryMd <- "<details style=\"margin: 0, 0, 0, 0; padding: 0, 0, 0, 0;\"><summary style=\"list-style-position: outside;\">"
        if (insertion[1] == 1) {
            countStr <- as.character(algTibble$count[i])
            countLen <- nchar(countStr)
            queryMd <- sprintf("%s%s<span style=\"letter-spacing: -0.3em;\">%s</span>", queryMd, sprintf("%9s", substr(countStr, 1, countLen - 1)), substr(countStr, countLen, countLen))
        } else {
            queryMd <- sprintf("%s%10d", queryMd, algTibble$count[i])
        }
        start <- 1
        refStart <- 1
        for (j in seq_len(length(insertion))) {
            refEnd <- insertionCollapse[j] - 1
            mdHigh <- snpHighlight(substr(algTibble$refNoGap[i], refStart, refEnd - 1), substr(algTibble$queryLine[i], start, insertion[j] - 2))
            refStart <- refEnd + 1
            lastBase <- substr(algTibble$queryLine[i], insertion[j] - 1, insertion[j] - 1)
            if (lastBase != toupper(substr(algTibble$refNoGap[i], refEnd, refEnd)) && lastBase != "-") {
                lastBase <- sprintf("<span style=\"color: red;\">%s</span>", lastBase)
            }
            queryMd <- sprintf("%s%s<span style=\"letter-spacing: -0.3em;\">%s|</span>", queryMd, mdHigh, lastBase)
            start <- insertion[j] + attr(insertion, which = "match.length")[j]
        }
        mdHigh <- snpHighlight(substr(algTibble$refNoGap[i], refStart, nchar(algTibble$refNoGap[i])), substr(algTibble$queryLine[i], start, nchar(algTibble$queryLine[i])))
        queryMd <- sprintf("%s%s</summary>%s</details>",
            queryMd,
            mdHigh,
            paste(insertionLines, collapse = '<br>')
        )
        allMd[2 * i] <- queryMd
    }
    return(allMd)
}
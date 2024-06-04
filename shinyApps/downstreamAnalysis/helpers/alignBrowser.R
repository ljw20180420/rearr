readAlgfiles <- function(algfiles) {
    algLiness <- vector(mode = "list", length = nrow(algfiles))
    for (i in seq_len(nrow(algfiles))) {
        algLiness[[i]] <- readLines(gzfile(algfiles$datapath[i]))
    }
    return(unlist(algLiness))
}

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

getMarkdownFromAlign <- function(algLines, algStart, algEnd) {
    allMd <- rep("", length(algLines) / 3 * 2)
    for (i in seq(algStart, algEnd)) {
        fields <- str_split_1(algLines[3 * i - 2], "\t")
        count <- fields[2]
        refId <- fields[4]
        cut1 <- as.integer(fields[16])
        reference <- algLines[3 * i - 1] |> str_replace_all("-", "")
        ref1Len <- gregexpr("[acgtn]", reference)[[1]][2]
        cut2 <- as.integer(fields[17]) - ref1Len
        refMd <- paste0(
            sprintf("<span style=\"background-color: lightgrey;\">%10s%s", refId, substr(reference, 1, cut1 - 1)),
            sprintf('<span style=\"letter-spacing: -0.3em;\">%s|</span>', substr(reference, cut1, cut1)),
            sprintf('<span style=\"color: red;\">%s</span>', substr(reference, cut1 + 1, ref1Len)),
            sprintf('<span style=\"color: blue;\">%s</span>', substr(reference, ref1Len + 1, ref1Len + cut2 - 1)),
            sprintf('<span style=\"color: blue; letter-spacing: -0.3em;\">%s</span><span style=\"letter-spacing: -0.3em;\">|</span>', substr(reference, ref1Len + cut2, ref1Len + cut2)),
            sprintf("%s</span>", substr(reference, ref1Len + cut2 + 1, nchar(reference)))
        )
        allMd[2 * i - 1] <- refMd
        query <- algLines[3 * i]
        insertion <- algLines[3 * i - 1] |> gregexpr(pattern = '-+') |> _[[1]]
        if (insertion[1] == -1) {
            queryMd <- sprintf("<details><summary style=\"list-style-position: outside;\">%10s%s</summary></details>", count, query)
            allMd[2 * i] <- queryMd
            next
        }
        insertionCollapse <- insertion - c(0, head(cumsum(attr(insertion, which = "match.length")), -1))
        insertionLines <- arrangeInsertion(query, insertion, insertionCollapse)
        queryMd <- "<details style=\"margin: 0, 0, 0, 0; padding: 0, 0, 0, 0;\"><summary style=\"list-style-position: outside;\">"
        if (insertion[1] == 1) {
            queryMd <- sprintf("%s%s<span style=\"letter-spacing: -0.3em;\">%s</span>", queryMd, sprintf("%9s", substr(count, 1, length(count) - 1)), substr(count, length(count), length(count)))
        } else {
            queryMd <- sprintf("%s%10s", queryMd, count)
        }
        start <- 1
        refStart <- 1
        for (j in seq_len(length(insertion))) {
            refEnd <- insertionCollapse[j] - 1
            mdHigh <- snpHighlight(substr(reference, refStart, refEnd - 1), substr(query, start, insertion[j] - 2))
            refStart <- refEnd + 1
            lastBase <- substr(query, insertion[j] - 1, insertion[j] - 1)
            if (lastBase != toupper(substr(reference, refEnd, refEnd)) && lastBase != "-") {
                lastBase <- sprintf("<span style=\"color: red;\">%s</span>", lastBase)
            }
            queryMd <- sprintf("%s%s<span style=\"letter-spacing: -0.3em;\">%s|</span>", queryMd, mdHigh, lastBase)
            start <- insertion[j] + attr(insertion, which = "match.length")[j]
        }
        mdHigh <- snpHighlight(substr(reference, refStart, nchar(reference)), substr(query, start, nchar(query)))
        queryMd <- sprintf("%s%s</summary>%s</details>",
            queryMd,
            mdHigh,
            paste(insertionLines, collapse = '<br>')
        )
        allMd[2 * i] <- queryMd
    }
    return(allMd)
}
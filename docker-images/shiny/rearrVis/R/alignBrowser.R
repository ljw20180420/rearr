arrangeInsertion <- function(query, insertion, insertionCollapse) {
  insertionLines <- tagList()
  for (i in seq_along(insertionCollapse)) {
    insertionAdded <- FALSE
    insertionSeg <- substr(
      query,
      insertion[i],
      insertion[i] + attr(insertion, which = "match.length")[i] - 1
    )
    for (j in seq_along(insertionLines)) {
      if (nchar(insertionLines[[j]]) < insertionCollapse[i]) {
        insertionLines[[j]] <- sprintf(
          "%s%s%s ",
          insertionLines[[j]],
          strrep(
            " ",
            insertionCollapse[i] - nchar(insertionLines[[j]]) - 1
          ),
          insertionSeg
        )
        insertionAdded <- TRUE
        break
      }
    }
    if (!insertionAdded) {
      insertionLines <- append(
        insertionLines,
        tagList(
          sprintf(
            "%s%s",
            strrep(" ", insertionCollapse[i] - 1),
            insertionSeg
          )
        )
      )
    }
  }
  insertionLines <- insertionLines |>
    lapply(function(x) {
      HTML(
        gsub(" ", "&nbsp;", paste0(strrep("#", 10), x))
      )
    })
  return(insertionLines)
}

snpHighlight <- function(refSeg, querySeg) {
  refArray <- stringr::str_split_1(toupper(refSeg), "")
  queryArray <- stringr::str_split_1(toupper(querySeg), "")
  diffs <- (refArray != queryArray & queryArray != "-") |> diff()
  diffs <- c(0, diffs, 0)
  boundaries <- diffs |> as.logical() |> which()
  if (length(boundaries) == 0) {
    return(tagList(querySeg))
  }
  mdStr <- tagList()
  start <- 1
  for (i in seq_len(length(boundaries) / 2)) {
    mdStr <- append(
      mdStr,
      tagList(
        substr(querySeg, start, boundaries[2 * i - 1] - 1),
        tags$span(
          style = "color: red;",
          substr(querySeg, boundaries[2 * i - 1], boundaries[2 * i] - 1)
        )
      )
    )
    start <- boundaries[2 * i]
  }
  return(append(mdStr, tagList(substr(querySeg, start, nchar(querySeg)))))
}

getMarkdownFromAlign <- function(algTibble) {
  allMd <- rep("", nrow(algTibble) * 2)
  for (i in seq_len(nrow(algTibble))) {
    refMd <- tags$span(
      style = "background-color: lightgrey;",
      sprintf(
        "%010d%s",
        algTibble$refId[i],
        substr(algTibble$refNoGap[i], 1, algTibble$cut1[i] - 1)
      ),
      tags$span(
        style = "letter-spacing: -0.3em;",
        sprintf(
          "%s|",
          substr(algTibble$refNoGap[i], algTibble$cut1[i], algTibble$cut1[i])
        )
      ),
      tags$span(
        style = "color: red;",
        substr(
          algTibble$refNoGap[i],
          algTibble$cut1[i] + 1,
          algTibble$ref1Len[i]
        )
      ),
      tags$span(
        style = "color: blue;",
        substr(
          algTibble$refNoGap[i],
          algTibble$ref1Len[i] + 1,
          algTibble$ref1Len[i] + algTibble$cut2[i] - 1
        )
      ),
      tags$span(
        style = "color: blue; letter-spacing: -0.3em;",
        substr(
          algTibble$refNoGap[i],
          algTibble$ref1Len[i] + algTibble$cut2[i],
          algTibble$ref1Len[i] + algTibble$cut2[i]
        )
      ),
      tags$span(style = "letter-spacing: -0.3em;", "|"),
      substr(
        algTibble$refNoGap[i],
        algTibble$ref1Len[i] + algTibble$cut2[i] + 1,
        nchar(algTibble$refNoGap[i])
      )
    )
    allMd[2 * i - 1] <- gsub(
      "\n",
      "",
      htmltools::doRenderTags(refMd, indent = FALSE)
    )
    insertion <- algTibble$refLine[i] |> gregexpr(pattern = "-+") |> _[[1]]
    if (insertion[1] == -1) {
      queryMd <- tags$details(
        tags$summary(
          style = "list-style-position: outside;",
          sprintf("%010d%s", algTibble$count[i], algTibble$queryLine[i])
        )
      )
      allMd[2 * i] <- gsub(
        "\n",
        "",
        htmltools::doRenderTags(queryMd, indent = FALSE)
      )
      next
    }
    insertionCollapse <- insertion -
      c(0, head(cumsum(attr(insertion, which = "match.length")), -1))
    insertionLines <- arrangeInsertion(
      algTibble$queryLine[i],
      insertion,
      insertionCollapse
    )

    queryMd <- tagList()
    if (insertion[1] == 1) {
      countStr <- as.character(algTibble$count[i])
      countLen <- nchar(countStr)

      queryMd <- append(
        queryMd,
        tagList(
          sprintf("%09s", substr(countStr, 1, countLen - 1)),
          tags$span(
            style = "letter-spacing: -0.3em;",
            substr(countStr, countLen, countLen)
          )
        )
      )
    } else {
      queryMd <- append(
        queryMd,
        tagList(sprintf("%010d", algTibble$count[i]))
      )
    }

    start <- 1
    refStart <- 1
    for (j in seq_along(insertion)) {
      refEnd <- insertionCollapse[j] - 1
      mdHigh <- snpHighlight(
        substr(algTibble$refNoGap[i], refStart, refEnd - 1),
        substr(algTibble$queryLine[i], start, insertion[j] - 2)
      )
      refStart <- refEnd + 1

      lastBase <- substr(
        algTibble$queryLine[i],
        insertion[j] - 1,
        insertion[j] - 1
      )

      queryMd <- append(queryMd, mdHigh)
      queryMd <- append(
        queryMd,
        tagList(
          tags$span(
            style = "letter-spacing: -0.3em;",
            ifelse(
              lastBase !=
                toupper(substr(algTibble$refNoGap[i], refEnd, refEnd)) &&
                lastBase != "-",
              tags$span(style = "color: red;", lastBase),
              lastBase
            ),
            "|",
          )
        )
      )

      start <- insertion[j] + attr(insertion, which = "match.length")[j]
    }
    mdHigh <- snpHighlight(
      substr(algTibble$refNoGap[i], refStart, nchar(algTibble$refNoGap[i])),
      substr(algTibble$queryLine[i], start, nchar(algTibble$queryLine[i]))
    )
    queryMd <- append(queryMd, mdHigh)

    queryMd <- do.call(
      tags$details,
      c(
        tagList(
          do.call(
            tags$summary,
            c(list(style = "list-style-position: outside;"), queryMd)
          )
        ),
        insertionLines
      )
    )

    allMd[2 * i] <- gsub(
      "\n",
      "",
      htmltools::doRenderTags(queryMd, indent = FALSE)
    )
  }
  return(allMd)
}

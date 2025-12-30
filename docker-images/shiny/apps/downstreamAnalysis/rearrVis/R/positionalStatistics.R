#' Right pad sequences
#'
#' Right pad sequences to the same length as the longest sequence by "-".
#'
#' @param seqs sequence vector
#' @return matrix with each row a right-padded sequence
#' @export
#'
#' @examples
#' refGapMat <- extendToSameLength(seqs)
extendToSameLength <- function(seqs) {
  maxLen <- max(nchar(seqs))
  return(
    seqs |>
      stringr::str_pad(width = maxLen, side = "right", pad = "-") |>
      strsplit("") |>
      unlist() |>
      matrix(ncol = maxLen, byrow = TRUE)
  )
}

#' Align sequence by cut site and pad on both sides
#'
#' All sequences are aligned according to the cut site. Then each sequence is padded by "-" on both sides such that the results sequences are of the same length.
#'
#' @param seqs sequence vector
#' @param cuts cut site vector
#' @return matrix with each row a padded sequence
#' @export
#'
#' @examples
#' refMat <- extendToAlignCut(seqs)
extendToAlignCut <- function(seqs, cuts) {
  leftMax <- max(cuts)
  rightMax <- max(nchar(seqs) - cuts)
  return(
    paste0(
      substr(seqs, 1, cuts) |>
        stringr::str_pad(width = leftMax, side = "left", pad = "-"),
      substr(seqs, cuts + 1, nchar(seqs)) |>
        stringr::str_pad(width = rightMax, side = "right", pad = "-")
    ) |>
      strsplit("") |>
      unlist() |>
      matrix(ncol = leftMax + rightMax, byrow = TRUE)
  )
}

#' Split character vector according to lengths
#'
#' Split a character vector into consecutive subvectors with specified length, and transform character subvectors to substrings.
#'
#' @param vec a character vector
#' @param strLens length of each substring
#' @return a vector of substrings
#' @export
#'
#' @examples
#' vectorToStringVector(c("A", "B", "C", "D", "E"), c(2, 3))
#' [1] "AB"  "CDE"
vectorToStringVector <- function(vec, strLens) {
  split(vec, rep.int(seq_along(strLens), strLens)) |>
    vapply(function(x) paste0(x, collapse = ""), "", USE.NAMES = FALSE)
}

#' Melt matrix
#'
#' Melt the matrix with position relative to cut site recorded in column pos and base type (rownames) recorded in column type. matrix values are recorded in column count.
#'
#' @param mat the base frequency matrix to melt
#' @param cut the cut site
#' @return melted matrix as tibble with columns count, pos and type
#' @export
#'
#' @examples
#' baseTibble <- posMatrixToTibble(baseFreq, maxCut)
posMatrixToTibble <- function(mat, cut) {
  tibble::tibble(
    count = c(mat),
    pos = rep(seq_len(ncol(mat)) - 0.5 - cut, each = nrow(mat)),
    type = factor(rep(rownames(mat), times = ncol(mat)), levels = rownames(mat))
  )
}

#' Stack bar plot of positional base frequencies or deletion/snp/match frequencies
#'
#' Base frequencies (including deletion frequencies expressed in "-") or deletion/snp/match frequencies at each position relative to cut site are visualized as stacked bar plot. Insertions are expressed by step plot with step height as the insertion frequencies and step position centered at where the insertions actually happens (between two base positions).
#'
#' @param inputTibble tibble with columns count (base frequencies), pos (relative position to cut site) and type (base type)
#' @param insertCount tibble with columns count (insertion frequencies) and pos (relative position to cut site)
#' @param refTempFile temp pdf file to save ggplot output
#' @return html iframe tag to be rendered by shiny::renderUI
#' @export
#'
#' @examples
#' iframeTag <- drawPositionalStatic(baseTibble, insertCount, tempFile)
#' iframeTag <- drawPositionalStatic(MSDTibble, insertCount, tempFile)
drawPositionalStatic <- function(inputTibble, insertCount, refTempFile) {
  ggFig <- inputTibble |>
    ggplot2::ggplot(ggplot2::aes(pos, count)) +
    ggplot2::geom_col(ggplot2::aes(fill = type), width = 1) +
    ggplot2::geom_step(
      ggplot2::aes(pos, count, color = "black"),
      data = insertCount[2:(nrow(insertCount) - 1), ],
      direction = "mid"
    ) +
    ggplot2::scale_x_continuous(
      name = "position relative to cut",
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_color_identity(
      name = NULL,
      guide = ggplot2::guide_legend(),
      labels = "insertion"
    )
  ggplot2::ggsave(
    refTempFile,
    plot = ggFig,
    height = 1200,
    width = 3600,
    unit = "px"
  )
  htmltools::tags$iframe(
    src = sub("^www/", "", refTempFile),
    height = "600px",
    width = "100%"
  )
}

#' Get query base frequencies
#'
#' The counts are grouped and summed with each column of queryMat as labels. The insertions (relative to reference line) in query lines are collapsed. Then the query lines are aligned according to cut site.
#'
#' @param queryMat a matrix of query lines
#' @param counts the counts for lines in queryMat
#' @return a matrix with five rows named "-", "A", "C", "G", "T". Each row records the positional frequency of each base.
#' @export
#'
#' @examples
#' baseFreq <- getPositionalBaseFreq(queryMat, counts)
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

#' Get query mutation frequencies
#'
#' The counts are grouped and summed with each column of queryMat as labels. The insertions (relative to reference line) in query lines are collapsed. Then the query lines are aligned according to cut site.
#'
#' @param queryMat a matrix of query lines
#' @param refMat a matrix of reference lines
#' @param counts the counts for lines in queryMat
#' @param cut the cut site
#' @return a tibble with columns count (frequencies of mutation), pos (relative position to cut site) and type (one of "delete", "SNP", "match").
#' @export
#'
#' @examples
#' MSDTibble <- getPositionalMSDTibble(queryMat, refMat, counts, maxCut)
getPositionalMSDTibble <- function(queryMat, refMat, counts, cut) {
  refMat <- toupper(refMat)
  MSDFreq <- rbind(
    colSums((queryMat == "-") * counts),
    colSums((queryMat != refMat & queryMat != "-") * counts),
    colSums((queryMat == refMat) * counts)
  )
  rownames(MSDFreq) <- c("delete", "SNP", "match")
  MSDFreq |> posMatrixToTibble(cut)
}

#' Get insertion frequencies
#'
#' Get positional insertion frequencies relative to the cut site.
#'
#' @param refList gapped reference lines in the two-line alignments
#' @param cuts cut sites
#' @param maxCutDown the maximal length of reference part downstream to the cut site
#' @return tibble with columns count (insertion frequencies) and pos (relative position to cut site)
#' @export
#'
#' @examples
#' insertCount <- calInsertionCount(refList, cuts, maxCutDown)
calInsertionCount <- function(refList, cuts, maxCutDown) {
  maxCut <- max(cuts)
  insertList <- vector("list", length(refList))
  for (i in seq_along(refList)) {
    mask <- refList[[i]] != "-"
    insertList[[i]] <- cumsum(mask)[!mask] - cuts[i]
  }
  histCount <- insertList |> unlist() |> table()
  fullCount <- rep(0, maxCut + maxCutDown + 1)
  fullCount[as.integer(names(histCount)) + maxCut + 1] <- histCount
  tibble::tibble(
    count = fullCount,
    pos = seq(-maxCut, maxCutDown)
  )
}

#' Adapter matrix to display resolution
#'
#' Because display resolution is only 1080 in height, matrix with large number of rows cannot be displayed smoothly and correctly. This function transforms matrix with large row numbers to a display matrix with specified row number (resolution). Each row of the matrix may repeat multiple times (counts).
#'
#' @param mat the matrix to display
#' @param counts the repeats of matrix rows
#' @param resolution the output matrix row number to fit the display resolution
#' @return description of the return value
#'
#' @examples
#' resolutionMatrix <- resolutionMatrix(mat, counts, resolution)
resolutionMatrix <- function(mat, counts, resolution) {
  totalCount <- sum(counts)
  resolutionCount <- round(seq(1, totalCount + 1, length.out = resolution + 1))
  resolutionMatrix <- matrix(0, resolution, ncol(mat))
  cumCount <- 0
  j <- 1
  for (i in seq_along(counts)) {
    cumCount <- cumCount + counts[i]
    if (cumCount < resolutionCount[j + 1]) {
      resolutionMatrix[j, ] <- resolutionMatrix[j, ] + mat[i, ] * counts[i]
    } else {
      resolutionMatrix[j, ] <- resolutionMatrix[j, ] +
        mat[i, ] * (resolutionCount[j + 1] - 1 + counts[i] - cumCount)
      j <- j + 1
      while (resolutionCount[j + 1] <= cumCount) {
        resolutionMatrix[j, ] <- mat[i, ] *
          (resolutionCount[j + 1] - resolutionCount[j])
        j <- j + 1
      }
      resolutionMatrix[j, ] <- mat[i, ] * (cumCount - resolutionCount[j] + 1)
    }
  }
  return(resolutionMatrix)
}

#' Display values stored in matrix rows
#'
#' Matrix rows with repeats are adapted to specified resolution. Then the matrix is melted for ploting by ggplot.
#'
#' @param mat matrix to process
#' @param counts repeats of matrix rows
#' @param resolution display resolution to adapted to
#' @param cut cut site
#' @return tibble with columns x (resolution matrix column index), y (resolution matrix row index), value (resolution matrix value) and height (resolution matrix pixel heights for rows)
#' @export
#'
#' @examples
#' resolutionTibble <- getPositionalReads(mat, counts, resolution, cut)
getPositionalReads <- function(mat, counts, resolution, cut) {
  resolutionMat <- resolutionMatrix(mat, counts, resolution)
  resolutionCount <- round(seq(1, sum(counts) + 1, length.out = resolution + 1))
  vertCent <- (resolutionCount[-length(resolutionCount)] +
    resolutionCount[-1]) /
    2 -
    1
  resolutionMat |>
    tibble::as_tibble(
      .name_repair = "universal_quiet",
      rownames = "pos1"
    ) |>
    tidyr::pivot_longer(
      cols = tidyselect::starts_with("..."),
      names_to = "pos2"
    ) |>
    dplyr::mutate(
      pos1 = as.integer(pos1),
      pos2 = as.integer(sub(r"(^\.\.\.)", "", pos2))
    ) |>
    dplyr::mutate(
      x = pos2 - cut - 0.5,
      y = vertCent[pos1],
      value = value,
      height = resolutionCount[pos1 + 1] - resolutionCount[pos1]
    ) |>
    dplyr::select(x, y, value, height)
}

#' Draw resolution tibble
#'
#' Resolution tibble from matrix adapted to specified resolution are visualized to ggplot2.
#'
#' @param tibb resolution tibble to plot
#' @param maxCut the maximal length of reference part upstream to cut site
#' @param maxCutDown the maximal length of reference part downstream to cut site
#' @param refTempFile temp pdf file to save ggplot output
#' @return html iframe tag to be rendered by shiny::renderUI
#' @export
#'
#' @examples
#' iframeTag <- drawPositionalReads(tibb, maxCut, maxCutdown, refTempFile)
drawPositionalReads <- function(tibb, maxCut, maxCutdown, refTempFile) {
  ggFig <- ggplot2::ggplot(
    tibb,
    ggplot2::aes(x = x, y = y, fill = value, height = height)
  ) +
    ggplot2::geom_tile(width = 1) +
    ggplot2::scale_fill_gradient(low = "white", high = "red") +
    ggplot2::scale_x_continuous(
      name = "position relative to cut",
      limits = c(-maxCut, maxCutdown),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(name = "reads", expand = c(0, 0))
  ggplot2::ggsave(
    refTempFile,
    plot = ggFig,
    height = 1200,
    width = 3600,
    unit = "px"
  )
  htmltools::tags$iframe(
    src = sub("^www/", "", refTempFile),
    height = "600px",
    width = "100%"
  )
}

#' Use ggseqlogo to represent frequency matrix
#'
#' Base frequencies (including deletion frequencies expressed in "-") at each position relative to cut site are visualized by ggseqlogo.
#'
#' @param baseFreq frequency matrix of base
#' @param method ggseqlogo method ("bits" or "probability") influence the heights of the logo
#' @param namespace alphabet of logo
#' @param refTempFile temp pdf file to save ggplot output
#' @return html iframe tag to be rendered by shiny::renderUI
#' @export
#'
#' @examples
#' iframeTag <- drawPositionalLogo(baseFreq, "bits", "ACGT", refTempFile)
drawPositionalLogo <- function(
  baseFreq,
  method,
  namespace,
  refTempFile
) {
  ggFig <- ggplot2::ggplot() +
    ggseqlogo::geom_logo(
      data = baseFreq,
      method = method,
      namespace = namespace
    ) +
    ggseqlogo::theme_logo()
  # Modify scale in existing ggplot object without replacing the scale. https://stackoverflow.com/questions/63691661/modify-scale-in-existing-ggplot-object-without-replacing-the-scale
  ggFig$scales$scales[[1]]$breaks <- NULL

  ggplot2::ggsave(
    refTempFile,
    plot = ggFig,
    height = 1200,
    width = 3600,
    unit = "px"
  )
  htmltools::tags$iframe(
    src = sub("^www/", "", refTempFile),
    height = "600px",
    width = "100%"
  )
}

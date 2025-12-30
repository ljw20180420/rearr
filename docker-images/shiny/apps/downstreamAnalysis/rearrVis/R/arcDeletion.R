#' Get deletion info
#'
#' Get the range of deletions for each query.
#'
#' @param algTibble tibble with alignment info
#' @return tibble with deletion info
#' @export
#'
#' @examples
#' arcDelTibble <- getArcDelTibble(algTibble)
getArcDelTibble <- function(algTibble) {
  lposes <- algTibble$refLine |>
    strsplit("") |>
    lapply(function(refLine) {
      c(0, cumsum(refLine != "-"))
    })
  delRegs <- algTibble$queryLine |> gregexpr(pattern = "-+")
  arcDelTibble <- tibble::tibble(
    count = algTibble$count,
    delStart1 = vector("list", length(delRegs)),
    delEnd1 = vector("list", length(delRegs)),
    delStart2 = vector("list", length(delRegs)),
    delEnd2 = vector("list", length(delRegs)),
  )
  for (i in seq_along(delRegs)) {
    if (delRegs[[i]][1] == -1) {
      next
    }
    delStart <- lposes[[i]][delRegs[[i]]]
    delEnd <- delStart + attributes(delRegs[[i]])$match.length
    mask1 <- delStart < algTibble$ref1Len[i]
    mask2 <- delEnd > algTibble$ref1Len[i]
    arcDelTibble$delStart1[[i]] <- delStart[mask1] - algTibble$cut1[i]
    arcDelTibble$delEnd1[[i]] <- pmin(delEnd[mask1], algTibble$ref1Len[i]) -
      algTibble$cut1[i]
    arcDelTibble$delStart2[[i]] <- pmax(delStart[mask2], algTibble$ref1Len[i]) -
      algTibble$ref1Len[i] -
      algTibble$cut2[i]
    arcDelTibble$delEnd2[[i]] <- delEnd[mask2] -
      algTibble$ref1Len[i] -
      algTibble$cut2[i]
  }
  return(arcDelTibble)
}

#' Express deletions by arcs
#'
#' Infer the arcs segments from insertion info.
#'
#' @param arcDelTibble tibble containing deletion info
#' @param segNum segment number to smooth the arc
#' @return tibble with arc segments
#' @export
#'
#' @examples
#' arcDelSegTibble <- getArcSegment(arcDelTibble, segNum)
getArcSegment <- function(arcDelTibble, segNum) {
  arcDelSegTibble <- arcDelTibble |>
    dplyr::mutate(
      x0 = (delStart + delEnd) / 2,
      r = (delEnd - delStart) / 2,
      x = vector("list", nrow(arcDelTibble)),
      y = vector("list", nrow(arcDelTibble)),
      xend = vector("list", nrow(arcDelTibble)),
      yend = vector("list", nrow(arcDelTibble)),
      count = count
    )

  for (i in seq_len(nrow(arcDelSegTibble))) {
    arcDelSegTibble$x[[i]] <- arcDelSegTibble$x0[i] +
      arcDelSegTibble$r[i] *
        cos(seq(0, pi * (segNum - 1) / segNum, length.out = segNum - 1))
    arcDelSegTibble$y[[i]] <- arcDelSegTibble$r[i] *
      sin(seq(0, pi * (segNum - 1) / segNum, length.out = segNum - 1))
    arcDelSegTibble$xend[[i]] <- arcDelSegTibble$x0[i] +
      arcDelSegTibble$r[i] * cos(seq(pi / segNum, pi, length.out = segNum - 1))
    arcDelSegTibble$yend[[i]] = arcDelSegTibble$r[i] *
      sin(seq(pi / segNum, pi, length.out = segNum - 1))
  }

  arcDelSegTibble |>
    select(x, y, xend, yend, count) |>
    unnest(c(x, y, xend, yend))
}

#' Plot deletion arcs
#'
#' Visualize deletions by arcs. The arc line width is the deletion count. The arc range is the deletion range.
#'
#' @param arcDelSegTibble tibble with arc segments
#' @param limits the limits of the reference coordinate
#' @param arcDeleteTempFile the pdf temp file
#' @return the html iframe tag to be displayed in shiny::renderUI
#' @export
#'
#' @examples
#' iframeTag <- plotArcDelTibble(arcDelSegTibble, limits, arcDeleteTempFile)
plotArcDelTibble <- function(arcDelSegTibble, limits, arcDeleteTempFile) {
  ggFig <- arcDelSegTibble |>
    ggplot2::ggplot(ggplot2::aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      linewidth = count
    )) +
    ggplot2::geom_segment(alpha = 0.1) +
    ggplot2::scale_linewidth_continuous(range = c(0.1, 2)) +
    ggplot2::scale_x_continuous(
      limits = limits,
      name = "pos",
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(name = NULL, expand = c(0, 0))
  ggplot2::ggsave(
    arcDeleteTempFile,
    plot = ggFig,
    height = 1200,
    width = 3600,
    unit = "px"
  )
  htmltools::tags$iframe(
    src = sub("^www/", "", arcDeleteTempFile),
    height = "600px",
    width = "100%"
  )
}

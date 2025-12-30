#' Get range tibble of micro-homology
#'
#' Record the ramge pf micro-homology in rows of tibble.
#'
#' @param ref1 the reference upstream to the cut site
#' @param ref2 the reference downstream to the cut site
#' @param cut1 the cut site in ref1
#' @param cut2 the cut site in ref2
#' @return a tibble with the four columns: pos1low (the row start position of mh), pos1up (the row end position of mh), shift (pos1up - pos2up), cls (the index of mh)
#' @export
#'
#' @examples
#' mhTibble <- getMicroHomologyTibble(ref1, ref2, cut1, cut2)
getMicroHomologyTibble <- function(ref1, ref2, cut1, cut2) {
  ref1vec <- stringr::str_split(toupper(ref1), "")[[1]]
  ref2vec <- stringr::str_split(toupper(ref2), "")[[1]]
  mh_matrix <- matrix(
    as.integer(
      rep(ref1vec, time = length(ref2vec)) ==
        rep(ref2vec, each = length(ref1vec))
    ),
    nrow = length(ref1vec)
  )
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
    mh_matrix |>
      tibble::as_tibble(.name_repair = "universal_quiet", rownames = "pos1") |>
      tidyr::pivot_longer(
        cols = tidyselect::starts_with("..."),
        names_to = "pos2",
        values_to = "cls"
      ) |>
      dplyr::mutate(
        pos1 = as.integer(pos1),
        pos2 = as.integer(sub(r"(^\.\.\.)", "", pos2))
      ) |>
      dplyr::filter(cls > 0) |>
      dplyr::summarise(
        pos1low = min(pos1) - 1 - cut1,
        pos1up = max(pos1) - cut1,
        shift = pos1up - max(pos2) + cut2,
        .by = "cls"
      )
  )
}

#' Visualize micro-homology
#'
#' Visualize the actual counts of editing output by heatmap and micro-homologies by paths passing the whole mhs.
#'
#' @param mhTibbleSub the tibble recording micro-homology ranges in rows
#' @param refEnd1Start2TibbleMicro the tibble containing actual counts of editing output for a certain reference
#' @param maxCut1 the maximal length of the ref1 part upstream to the cut site
#' @param maxCut2 the maximal length of the ref2 part upstream to the cut site
#' @param maxCut1down the maximal length of the ref1 part downstream to the cut site
#' @param maxCut2down the maximal length of the ref2 part downstream to the cut site
#' @param mode if set to "separate", then meaningly distribute the actual counts of editing output to the all micro-homology equivalence
#' @param mhMatrixTempFile the pdf tempfile to save ggplot result
#' @return the html iframe tag to be displayed by shiny::renderUI
#' @export
#'
#' @examples
#' iframeTag <- drawMicroHomologyHeatmap(mhTibbleSub, refEnd1Start2TibbleMicro, maxCut1, maxCut2, maxCut1down, maxCut2down, mode, mhMatrixTempFile)
drawMicroHomologyHeatmap <- function(
  mhTibbleSub,
  refEnd1Start2TibbleMicro,
  maxCut1,
  maxCut2,
  maxCut1down,
  maxCut2down,
  mode,
  mhMatrixTempFile
) {
  if (mode == "separate") {
    refEnd1Start2TibbleMicro <- dplyr::bind_rows(
      refEnd1Start2TibbleMicro |> dplyr::filter(cls == 0),
      refEnd1Start2TibbleMicro |>
        dplyr::filter(cls != 0) |>
        dplyr::group_by(cls) |>
        dplyr::mutate(count = count / dplyr::n()) |>
        dplyr::ungroup()
    )
  }
  mhPosTibble <- mhTibbleSub |>
    dplyr::mutate(pos2low = pos1low - shift, pos2up = pos1up - shift) |>
    dplyr::select(pos1low, pos1up, pos2low, pos2up)
  ggFig <- ggplot2::ggplot(refEnd1Start2TibbleMicro) +
    ggplot2::geom_tile(
      ggplot2::aes(x = pos2, y = pos1, fill = log10(count + 1)),
      height = 1,
      width = 1
    ) +
    # ggplot2::geom_segment(
    #   ggplot2::aes(
    #     x = pos2low,
    #     y = pos1low,
    #     xend = pos2up,
    #     yend = pos1up
    #   ),
    #   data = mhPosTibble
    # ) +
    ggplot2::scale_x_continuous(
      limits = c(-maxCut2 - 1, maxCut2down + 1),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-maxCut1 - 1, maxCut1down + 1),
      expand = c(0, 0)
    ) +
    ggplot2::scale_fill_gradientn(
      limits = c(0, NA),
      colors = c("blue", "white", "red")
    ) +
    ggplot2::scale_size_area(max_size = 2) +
    ggplot2::coord_equal(ratio = 1)
  ggplot2::ggsave(
    mhMatrixTempFile,
    plot = ggFig,
    height = 3600,
    width = 3600,
    unit = "px"
  )
  htmltools::tags$iframe(
    src = sub("^www/", "", mhMatrixTempFile),
    height = "1200px",
    width = "100%"
  )
}

#' Get editing output counts
#'
#' Get the actual counts of editing output for a specific reference.
#'
#' @param algTibble the tibble with columns ref1End (the reference end position of the first alignment block), cut1 (the cut site in ref1), ref2Start (the reference start position of the second alignment block), cut2 (the cut site in ref2), refId (the reference id), count (the count of the query)
#' @param microRefId the reference id to harvest
#' @return the tibble with the columns pos1 (the end position of the upstream alignment block relative to the cut site), pos2 (the end position of the downstream alignment block relative to the cut site), count (the count of the query), shift (pos1 - pos2)
#' @export
#'
#' @examples
#' refEnd1Start2Tibble <- getRefEnd1Start2Tibble(algTibble, microRefId)
getRefEnd1Start2Tibble <- function(algTibble, microRefId) {
  algTibble |>
    dplyr::mutate(
      pos1 = ref1End - cut1,
      pos2 = ref2Start - cut2,
      refId = refId,
      count = count,
    ) |>
    dplyr::filter(refId == microRefId) |>
    dplyr::summarise(count = sum(count), .by = c("pos1", "pos2")) |>
    dplyr::mutate(shift = pos1 - pos2)
}

#' Duplicate the editing output count upon micro-homology equivalence
#'
#' For each editing output, if it is located in a micro-homology, then duplicate its count to all equivalence editing output in the micro-homology range.
#'
#' @param refEnd1Start2Tibble the tibble recording the counts of editing outputs
#' @param mhTibbleSub the tibble recording micro-homology ranges in rows
#' @return the tibble duplicates the counts of editing outputs upon micro-homology equivalence
#' @export
#'
#' @examples
#' refEnd1Start2TibbleMicro <- getRefEnd1Start2TibbleMicro(refEnd1Start2Tibble, mhTibbleSub)
getRefEnd1Start2TibbleMicro <- function(refEnd1Start2Tibble, mhTibbleSub) {
  joinTibble <- refEnd1Start2Tibble |>
    dplyr::left_join(mhTibbleSub, by = "shift", relationship = "many-to-many")
  outRangeTibble <- joinTibble |>
    dplyr::summarise(
      inRange = any(pos1 >= pos1low & pos1 <= pos1up),
      count = dplyr::first(count),
      .by = c("pos1", "pos2")
    ) |>
    dplyr::filter(is.na(inRange) | !inRange) |>
    dplyr::mutate(inRange = NULL, cls = 0)
  inRangeTibble <- joinTibble |> dplyr::filter(pos1 >= pos1low, pos1 <= pos1up)
  if (nrow(inRangeTibble) == 0) {
    return(outRangeTibble)
  }
  inRangeTibble |>
    dplyr::rowwise() |>
    dplyr::mutate(pos1 = list(seq(pos1low, pos1up))) |>
    dplyr::ungroup() |>
    tidyr::unnest(pos1) |>
    dplyr::mutate(pos2 = pos1 - shift) |>
    dplyr::select(pos1, pos2, count, cls) |>
    dplyr::bind_rows(outRangeTibble)
}

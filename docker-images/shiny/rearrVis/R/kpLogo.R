getKpLogoAlgTarget <- function(
  algTibble,
  sgRNAs,
  editTarget,
  countThres,
  method
) {
  algTarget <- algTibble |>
    dplyr::mutate(target = editTarget) |>
    dplyr::summarise(
      count = sum(count),
      targetCount = sum(count * target),
      .by = refId
    ) |>
    dplyr::filter(count > countThres) |>
    dplyr::mutate(sgRNA = sgRNAs[refId + 1]) |>
    dplyr::select(-"refId")
  if (method == "weight") {
    algTarget <- algTarget |>
      dplyr::mutate(sgRNA = sgRNA, weight = targetCount / count) |>
      dplyr::select(sgRNA, weight)
  }
  return(algTarget)
}

plotKpLogoAlgTarget <- function(
  algTarget,
  method,
  region,
  kmer,
  outputKpLogoTempFile,
  weightKpLogoTempFile,
  targetKpLogoTempFile,
  bgFileKpLogoTempFile
) {
  if (method == "weight") {
    algTarget |> readr::write_tsv(weightKpLogoTempFile, col_names = FALSE)
    system2(
      "kpLogo",
      args = c(
        weightKpLogoTempFile,
        "-o",
        outputKpLogoTempFile,
        "-region",
        paste(region[1], region[2], sep = ","),
        "-weighted",
        "-k",
        kmer
      )
    )
  } else if (method == "background") {
    algTarget |>
      dplyr::select(sgRNA, targetCount) |>
      tidyr::uncount(targetCount) |>
      readr::write_tsv(targetKpLogoTempFile, col_names = FALSE)
    algTarget |>
      dplyr::select(sgRNA, count) |>
      tidyr::uncount(count) |>
      readr::write_tsv(bgFileKpLogoTempFile, col_names = FALSE)
    system2(
      "kpLogo",
      args = c(
        targetKpLogoTempFile,
        "-o",
        outputKpLogoTempFile,
        "-region",
        paste(region[1], region[2], sep = ","),
        "-bgfile",
        bgFileKpLogoTempFile,
        "-k",
        kmer
      )
    )
  }
  return(htmltools::tags$iframe(
    src = paste0(sub("^www/", "", outputKpLogoTempFile), ".all.pdf"),
    height = "1000px",
    width = "100%"
  ))
}

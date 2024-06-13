getKpLogoAlgTarget <- function(algTibble, sgRNAs, editTarget, countThres, method) {
    algTarget <- algTibble |> mutate(target = editTarget) |> summarise(count = sum(count), targetCount = sum(count * target), .by = refId) |> filter(count > countThres) |> mutate(sgRNA = sgRNAs[refId + 1]) |> select(-"refId")
    if (method == "weight") {
        algTarget <- algTarget |> mutate(sgRNA = sgRNA, weight = targetCount/count) |> select(sgRNA, weight)
    }
    return(algTarget)
}

plotKpLogoAlgTarget <- function(algTarget, method, region, kmer, outputKpLogo, weightKpLogo, targetKpLogo, bgFileKpLogo) {
    if (method == "weight") {
        algTarget |> write_tsv(weightKpLogo, col_names = FALSE)
        system2("kpLogo", args = c(
            weightKpLogo, "-o", outputKpLogo, "-region", paste(region[1], region[2], sep = ","), "-weighted", "-k", kmer
        ))
    } else if (method == "background") {
        algTarget |> select(sgRNA, targetCount) |> uncount(targetCount) |> write_tsv(targetKpLogo, col_names = FALSE)
        algTarget |> select(sgRNA, count) |> uncount(count) |> write_tsv(bgFileKpLogo, col_names = FALSE)
        system2("kpLogo", args = c(
            targetKpLogo, "-o", outputKpLogo, "-region", paste(region[1], region[2], sep = ","), "-bgfile", bgFileKpLogo, "-k", kmer
        ))
    }
    return(tags$iframe(src = paste0(sub("^www", "", outputKpLogo), ".all.pdf"), height = "1000px", width = "100%"))
}
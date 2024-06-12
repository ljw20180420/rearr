getArcDelTibble <- function(refLines, queryLines, counts, ref1Lens, cuts1, cuts2) {
    lposes <- refLines |> strsplit("") |> lapply(function(refLine){c(0, cumsum(refLine != "-"))})
    delRegs <- queryLines |> gregexpr(pattern = "-+")
    arcDelTibble <- tibble(
        count = counts,
        delStart1 = vector("list", length(delRegs)),
        delEnd1 = vector("list", length(delRegs)),
        delStart2 = vector("list", length(delRegs)),
        delEnd2 = vector("list", length(delRegs)),
    )
    for (i in seq_len(length(delRegs))) {
        if (delRegs[[i]][1] == -1) {
            next
        }
        delStart <- lposes[[i]][delRegs[[i]]]
        delEnd <- delStart + attributes(delRegs[[i]])$match.length
        mask1 <- delStart < ref1Lens[i]
        mask2 <- delEnd > ref1Lens[i]
        arcDelTibble$delStart1[[i]] <- delStart[mask1] - cuts1[i]
        arcDelTibble$delEnd1[[i]] <- pmin(delEnd[mask1], ref1Lens[i]) - cuts1[i]
        arcDelTibble$delStart2[[i]] <- pmax(delStart[mask2], ref1Lens[i]) - ref1Lens[i] - cuts2[i]
        arcDelTibble$delEnd2[[i]] <- delEnd[mask2] - ref1Lens[i] - cuts2[i]
    }
    return(arcDelTibble)
}

plotArcDelTibble <- function(arcDelTibble, limits) {
    arcDelTibble |>
        ggplot(aes(x0 = (delStart + delEnd) / 2, y0 = 0, r = (delEnd - delStart) / 2, start = - pi / 2, end = pi / 2)) +
        geom_arc(aes(linewidth = count), alpha = 0.1) +
        scale_linewidth_continuous(range = c(0.1, 2)) +
        scale_x_continuous(limits = limits, name = "pos", expand = c(0, 0)) +
        scale_y_continuous(name = NULL, expand = c(0, 0))
}
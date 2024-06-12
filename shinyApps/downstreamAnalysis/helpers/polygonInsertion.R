getPolyInsTibble <- function(refLines, counts, ref1Lens, cuts1, cuts2) {
    insRegs <- gregexpr("-+", refLines)
    insLens <- lapply(
        insRegs,
        function(insReg) {
            if (insReg[1] == -1) {
                return(NULL)
            }
            attributes(insReg)$match.length
        }
    )
    polyInsTibble <- tibble(
        count = counts,
        insPos1 = vector("list", length(insRegs)),
        insPos2 = vector("list", length(insRegs)),
        insLen1 = vector("list", length(insRegs)),
        insLen2 = vector("list", length(insRegs))
    )
    for (i in seq_len(length(insRegs))) {
        if (insRegs[[i]][1] == -1) {
            next
        }
        clen <- cumsum(attributes(insRegs[[i]])$match.length)
        insPos <- insRegs[[i]] - 1 - c(0, clen[seq_len(length(clen) - 1)])
        insLen <- attributes(insRegs[[i]])$match.length
        mask1 <- insPos <= ref1Lens[i]
        mask2 <- insPos >= ref1Lens[i]
        polyInsTibble$insPos1[[i]] <- insPos[mask1] - cuts1[i]
        polyInsTibble$insPos2[[i]] <- insPos[mask2] - ref1Lens[i] - cuts2[i]
        polyInsTibble$insLen1[[i]] <- insLen[mask1]
        polyInsTibble$insLen2[[i]] <- insLen[mask2]
    }
    return(polyInsTibble)
}

getPolyXY <- function(polyInsTibble, mode) {
    polyX <- rep(polyInsTibble$insPos, each = 4)
    if (mode == "down") {
        polyX[seq(3, length(polyX), 4)] <- polyX[seq(3, length(polyX), 4)] + polyInsTibble$insLen
    } else if (mode == "up") {
        polyX[seq(3, length(polyX), 4)] <- polyX[seq(3, length(polyX), 4)] - polyInsTibble$insLen
    } else {
        stop("mode must be down or up")
    }
    polyY <- rep(0, time = nrow(polyInsTibble) * 4)
    polyY[2:3 + rep(seq(0, length(polyY) - 1, 4), each = 2)] <- rep(polyInsTibble$count, each = 2)
    return(tibble(x = polyX, y = polyY))
}

plotPolyInsTibble <- function(polyXY, limits) {
    polyXY |>
        ggplot(aes(x, y)) +
        geom_polygon(color = "black", fill = NA, linewidth = 0.1) +
        scale_x_continuous(limits = limits, name = "pos", expand = c(0, 0)) +
        scale_y_continuous(name = "count", expand = c(0, 0))
}
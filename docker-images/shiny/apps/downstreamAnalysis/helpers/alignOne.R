alignOne <- function(ref1, ref2, cut1, cut2, correct, query, refFile, correctFile, algFile) {
    writeLines(
        sprintf(
            "0\t%s\t%d\t%d\t%s\t%d\n",
            ref1, cut1, cut2, ref2, nchar(ref2)
        ),
        con = refFile
    )
    writeLines(
        sprintf("%s\n", correct),
        con = correctFile
    )
    alignPipe = pipe(sprintf(
            'rearrangement 3<%s -u -3 -v -9 -s0 -6 -s1 4 -s2 2 -qv -5 | gawk -f correct_micro_homology.awk -- %s %s | tail -n+2 >%s',
            refFile,
            refFile,
            correctFile,
            algFile
    ))
    sprintf("%s\t1\t0", query) |> writeLines(con = alignPipe)
    return(readLines(algFile) |> paste(collapse = "<br>"))
}
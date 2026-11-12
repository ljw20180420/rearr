#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

strComp <- function(X) {
  return(seqinr::c2s(rev(seqinr::comp(seqinr::s2c(X)))))
}

ref_fasta <- args[1]
ref_out <- args[2]
gR <- args[3]

ref <- seqinr::read.fasta(ref_fasta, as.string = TRUE)

rvComp_seq <- strComp(ref[[1]][1])

align <- pwalign::pairwiseAlignment(
  toupper(gR),
  toupper(ref[[1]][1]),
  type = "global-local"
)
alignRevComp <- pwalign::pairwiseAlignment(
  toupper(gR),
  toupper(rvComp_seq),
  type = "global-local"
)

if (BiocGenerics::score(align) >= BiocGenerics::score(alignRevComp)) {
  seqinr::write.fasta(
    names = names(ref),
    sequences = ref[[1]][1],
    file.out = ref_out
  )
} else {
  seqinr::write.fasta(
    names = names(ref),
    sequences = rvComp_seq,
    file.out = ref_out
  )
}

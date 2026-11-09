#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
sgRNA <- Biostrings::DNAString(args[1])
bam <- normalizePath(args[2])
sgstart <- strtoi(args[3])

gdl <- GenomicRanges::GRanges(
  seqnames = S4Vectors::Rle(values = "ref", lengths = 1),
  ranges = IRanges::IRanges(start = sgstart + 1, end = sgstart + length(sgRNA)),
  name = S4Vectors::Rle(values = "ref", lengths = 1),
  strand = S4Vectors::Rle(values = BiocGenerics::strand("+"), lengths = 1),
  score = 0
)
crispr_set <- CrispRVariants::readsToTarget(
  reads = bam,
  target = gdl,
  reference = sgRNA,
  target.loc = 17
)
alns <- crispr_set$crispr_runs[[1]]$alns
for (i in seq_along(alns)) {
  cat(sprintf(
    fmt = "%s\t%s\n",
    names(alns)[i],
    GenomicAlignments::cigar(alns)[i]
  ))
}

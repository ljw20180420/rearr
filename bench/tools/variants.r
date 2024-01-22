#!/usr/bin/env Rscript
library("CrispRVariants")
library("Rsamtools")
library("GenomicAlignments")

args <- commandArgs(trailingOnly = TRUE)
ref <- Biostrings::DNAString(args[1])
bam <- args[2]
sgstart <- strtoi(args[3])
sglen <- strtoi(args[4])

gdl <- GRanges(seqnames = Rle("ref", 1), ranges = IRanges(sgstart + 1, sgstart + sglen), name = Rle("ref", 1), strand = Rle(strand("+"), 1), score = 0)
crispr_set <- readsToTarget(bam, gdl, reference = ref[(sgstart + 1):(sgstart + sglen)], target.loc = 17)
alns <- crispr_set$crispr_runs[[1]]$alns
for (i in seq_along(alns)){
  cat(sprintf("%s\t%s\n", names(alns)[i], cigar(alns)[i]))
}

#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config <- args[1]
fastq_folder <- args[2]
alignment_file <- args[3]

aln <- amplican::amplicanAlign(
  config = normalizePath(config),
  fastq_folder = normalizePath(fastq_folder),
  fastqfiles = 1
)

amplican::writeAlignments(
  aln,
  file = normalizePath(alignment_file),
  aln_format = "txt"
)

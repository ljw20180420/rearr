library(amplican)

args <- commandArgs(trailingOnly = TRUE)
config <- args[1]
fastq_folder <- args[2]
results_folder <- args[3]
amplicanPipeline(config, fastq_folder, results_folder, knit_reports = FALSE, fastqfiles = 1)

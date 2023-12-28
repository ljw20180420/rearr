library(amplican)

args <- commandArgs(trailingOnly = TRUE)
config <- args[1]
fastq_folder <- args[2]
results_folder <- args[3]

config <- normalizePath(config)
fastq_folder <- normalizePath(fastq_folder)
results_folder <- normalizePath(results_folder)
message("Checking write access...")
resultsFolder <- file.path(results_folder, "alignments")
if (!dir.exists(resultsFolder)) {
    dir.create(resultsFolder)
}
rds_file <- file.path(resultsFolder, "AlignmentsExperimentSet.rds")
aln <- amplicanAlign(config = config, fastq_folder = fastq_folder, fastqfiles = 1)
message("Saving alignments...")
writeAlignments(aln, file.path(resultsFolder, "alignments.txt"), "txt")

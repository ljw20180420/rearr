#!/usr/bin/env Rscript
# Usage
# rearr_render.r input.alg.cut1.cut2.ref1len

library(this.path)
args <- commandArgs(trailingOnly = TRUE)
algfile <- fs::path_abs(args[1])
rmarkdown::render(
    input = file.path(this.dir(), "tools/draw_figures.Rmd"),
    output_format = "html_document",
    output_file = paste(basename(algfile), "html", sep = "."),
    output_dir = dirname(algfile),
    intermediates_dir = dirname(algfile),
    knit_root_dir = dirname(algfile),
    params = list(algfile = algfile)
)
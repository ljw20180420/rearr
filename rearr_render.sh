#!/bin/bash
# Usage
# rearr_render.sh input.alg.cut1.cut2.ref1len

algfile="$(realpath $1)"
project_path="$(dirname $(realpath $0))"
Rscript --verbose -e '
    args = commandArgs(trailingOnly = TRUE)
    rmarkdown::render(
        input = file.path(args[1], "tools/draw_figures.Rmd"),
        output_format = "html_document",
        output_file = paste(basename(args[2]), "html", sep = "."),
        output_dir = dirname(args[2]),
        intermediates_dir = dirname(args[2]),
        knit_root_dir = dirname(args[2]),
        params = list(algfile = args[2])
    )
' $project_path $algfile
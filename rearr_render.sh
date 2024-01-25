#!/bin/bash
# Usage
# rearr_render.sh input.alg.cut1.cut2.ref1len

algfile="$(realpath $1)"
project_path="$(dirname $(realpath $0))"
quarto render "$project_path/tools/draw_figures.qmd" -P "algfile:$algfile"
mv "$project_path/tools/draw_figures.html" "$algfile.html"
rm -r "$(dirname $algfile)/draw_figures_files"
mv "$project_path/tools/draw_figures_files" "$(dirname $algfile)"
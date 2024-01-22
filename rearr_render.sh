#!/bin/bash
# Usage
# rearr_render.sh input.alg.cut1.cut2.ref1len

algfile="$(realpath $1)"
script_path=$(dirname $0)
quarto render "$script_path/tools/draw_figures.qmd" -P "algfile:$algfile"
mv "$script_path/tools/draw_figures.html" "$algfile.html"
rm -r "$(dirname $algfile)/draw_figures_files"
mv "$script_path/tools/draw_figures_files" "$(dirname $algfile)"
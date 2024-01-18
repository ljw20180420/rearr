#!/bin/bash
# Usage
# rearr_render.sh input.alg.cut1.cut2.ref1len

algfile="$(realpath $1)"
quarto render "$(dirname $0)/tools/draw_figures.qmd" -P "algfile:$algfile"
mv "$(dirname $0)/tools/draw_figures.html" "$algfile.html"
rm -r "$(dirname $algfile)/draw_figures_files"
mv "$(dirname $0)/tools/draw_figures_files" "$(dirname $algfile)"
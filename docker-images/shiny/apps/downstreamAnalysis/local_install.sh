#!/bin/bash

# Change directory to the script directory.
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

mkdir -p library
R -e 'withr::with_libpaths("./library", devtools::install("./rearrVis", dependencies = FALSE))'
R -e 'withr::with_libpaths("./library", devtools::install("./ggseqlogo", dependencies = FALSE))'

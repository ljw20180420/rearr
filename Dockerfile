# Analysis of CRISPR/Cas sequencing to get editing events
FROM ubuntu:22.04
LABEL maintainer="ljw2017@sjtu.edu.cn"
# install system enviroments
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends python3 python3-pip python3-dev r-base build-essential cmake bash gawk
# install python packages
RUN pip install -U pip && pip install biopython more-itertools numpy pysam scipy flask
# install R packages
RUN apt-get install -y --no-install-recommends r-cran-jsonlite r-cran-tidyverse r-cran-ggforce r-cran-patchwork r-cran-reticulate r-cran-ggtext r-cran-cairo r-cran-ggseqlogo && Rscript -e 'install.packages(c("this.path", "waffle"))'
# set working directory
WORKDIR /app
# copy aligner c++ source code
COPY ["Rearrangement/CMakeLists.txt", "Rearrangement/main.cpp", "/app/Rearrangement/"]
COPY ["Rearrangement/headers", "/app/Rearrangement/headers"]
# build aligner c++ source code
RUN cmake -DCMAKE_BUILD_TYPE=Release -S Rearrangement -B Rearrangement/build && make -C Rearrangement/build
# copy scripts
COPY ["./rearr_run.sh", "./rearr_render.sh", "/app/"]
COPY ["./tools/correct_micro_homology.AWK", "/app/tools/"]
COPY ["./tools/draw_figures_docker.Rmd", "/app/tools/draw_figures.Rmd"]

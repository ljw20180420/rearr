# Analysis of CRISPR/Cas sequencing to get editing events
FROM ubuntu:22.04
LABEL maintainer="ljw2017@sjtu.edu.cn"
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends python3 python3-pip r-base bash cmake build-essential gawk
RUN pip install -U pip && pip install biopython more-itertools numpy pysam scipy flask
RUN apt-get install -y --no-install-recommends r-cran-jsonlite r-cran-tidyverse r-cran-ggforce r-cran-patchwork r-cran-reticulate r-cran-ggtext r-cran-cairo r-cran-ggseqlogo
RUN Rscript -e 'install.packages(c("this.path", "waffle"))'
WORKDIR /app
COPY ["Rearrangement/CMakeLists.txt", "Rearrangement/main.cpp", "/app/Rearrangement/"]
COPY ["Rearrangement/headers", "/app/Rearrangement/headers"]
RUN cmake -DCMAKE_BUILD_TYPE=Release -S Rearrangement -B Rearrangement/build && make -C Rearrangement/build
COPY ["./rearr_run.sh", "./rearr_render.sh", "/app/"]
COPY ["./tools/correct_micro_homology.AWK", "./tools/draw_figures.Rmd", "/app/tools/"]
RUN apt-get install -y --no-install-recommends python3-dev

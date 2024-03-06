# TODO: not use root user
# TODO: use minimal executable

# Analysis of CRISPR/Cas sequencing to get editing events
FROM ubuntu:22.04
LABEL maintainer="ljw2017@sjtu.edu.cn"
# install system enviroments
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends python3 python3-pip python3-dev r-base build-essential libcurl4-openssl-dev libssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev cmake bash gawk bsdmainutils libncurses5-dev libncursesw5-dev pandoc
# install python packages
RUN pip install -U pip && pip install biopython more-itertools numpy pysam scipy flask
# install R packages
RUN apt-get install -y --no-install-recommends r-cran-jsonlite r-cran-ggforce r-cran-reticulate r-cran-ggtext r-cran-cairo r-cran-ggseqlogo && Rscript -e 'install.packages(c("tidyverse", "patchwork", "this.path", "waffle"))'
# set working directory
WORKDIR /app
# copy aligner c++ source code
COPY ["Rearrangement/CMakeLists.txt", "Rearrangement/main.cpp", "/app/Rearrangement/"]
COPY ["Rearrangement/headers", "/app/Rearrangement/headers"]
# build aligner c++ source code
RUN cmake -DCMAKE_BUILD_TYPE=Release -S Rearrangement -B Rearrangement/build && make -C Rearrangement/build
# build less
COPY ["./less-643.zip", "/app/"]
RUN unzip less-643.zip && cd "less-643" && ./configure && make
# build kpLogo
COPY ["kpLogo-1.1.zip", "/app/barcode/"]
RUN unzip barcode/kpLogo-1.1.zip -d barcode && make -C "barcode/kpLogo-1.1/src"
# build pv
COPY ["pv-1.8.5.zip", "/app/"]
RUN unzip pv-1.8.5.zip && cd pv-1.8.5 && ./configure && make
# copy scripts for single target alignment
COPY ["rearr_run.sh", "rearr_render.sh", "rearr_view.sh", "/app/"]
COPY ["tools/correct_micro_homology.AWK", "tools/align_align.py", "tools/draw_figures.Rmd", "/app/tools/"]
# copy barcode library

# TODO: not use root user
# TODO: use minimal executable

# build static binary tools for ubuntu:22.04
FROM ubuntu:22.04 AS ubuntu-bin
LABEL maintainer="ljw2017@sjtu.edu.cn"
## set working directory
WORKDIR /app
## install system enviroments for build
RUN apt-get update && apt-get install -y --no-install-recommends unzip build-essential libncurses5-dev
## build less
COPY ["./less-643.zip", "/app/"]
RUN unzip less-643.zip 
RUN cd "less-643" && ./configure LDFLAGS="-static" && make
# install python
# RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-pip
# RUN pip install -U pip && pip install pyinstaller numpy scipy flask
# # install system enviroments
# RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends python3-dev r-base libcurl4-openssl-dev libssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev cmake bash gawk bsdmainutils libncursesw5-dev pandoc samtools bedtools bowtie2 cutadapt ghostscript imagemagick
# # fix policy problem of imagemagick
# RUN sed -i -r 's#<policy domain="coder" rights="none" pattern="(PS|PS2|PS3|EPS|PDF|XPS)" />#<!-- policy domain="coder" rights="none" pattern="\1" />#g' /etc/ImageMagick-6/policy.xml
# # install R packages
# RUN apt-get install -y --no-install-recommends r-cran-jsonlite r-cran-ggforce r-cran-reticulate r-cran-ggtext r-cran-cairo r-cran-ggseqlogo && Rscript -e 'install.packages(c("tidyverse", "patchwork", "this.path", "waffle"))'
# # copy aligner c++ source code
# COPY ["Rearrangement/CMakeLists.txt", "Rearrangement/main.cpp", "/app/Rearrangement/"]
# COPY ["Rearrangement/headers", "/app/Rearrangement/headers"]
# # build aligner c++ source code
# RUN cmake -DCMAKE_BUILD_TYPE=Release -S Rearrangement -B Rearrangement/build && make -C Rearrangement/build
# # build kpLogo
# COPY ["kpLogo-1.1.zip", "/app/barcode/"]
# RUN unzip barcode/kpLogo-1.1.zip -d barcode && make -C "barcode/kpLogo-1.1/src"
# # build pv
# COPY ["pv-1.8.5.zip", "/app/"]
# RUN unzip pv-1.8.5.zip && cd pv-1.8.5 && ./configure && make
# # copy scripts for single target alignment
# COPY ["rearr_run.sh", "rearr_render.r", "rearr_view.sh", "/app/"]
# COPY ["tools/correct_micro_homology.AWK", "tools/align_align.py", "tools/draw_figures.Rmd", "/app/tools/"]
# # copy barcode library
# COPY ["barcode/sx/get*.sh", "barcode/sx/index_spliter.sh", "barcode/sx/infer_csvfile.sh", "/app/barcode/sx/"]
# COPY ["barcode/tools/demultiplex.sh", "barcode/tools/barcode_align.sh", "barcode/tools/barcode_post_process.sh", "barcode/tools/run_kpLogo.sh", "barcode/tools/trouble_shooting.r", "barcode/tools/barcode_render.sh", "/app/barcode/tools/"]

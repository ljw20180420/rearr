# Analysis of CRISPR/Cas sequencing to get editing events
FROM ubuntu:22.04
LABEL maintainer="ljw2017@sjtu.edu.cn"
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends python3 python3-pip r-base bash
RUN pip install -U pip && pip install biopython more-itertools numpy pysam scipy flask
RUN apt-get install -y --no-install-recommends r-cran-jsonlite r-cran-tidyverse r-cran-ggforce r-cran-patchwork r-cran-reticulate r-cran-ggtext r-cran-cairo r-cran-ggseqlogo
RUN Rscript -e 'install.packages(c("this.path", "waffle"))'
ENV APPROOT="/app"
WORKDIR $APPROOT
COPY ["./rearr_run.sh", "${APPROOT}"]

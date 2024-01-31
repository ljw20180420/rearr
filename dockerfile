FROM ubuntu:20.04

USER root
WORKDIR /root

SHELL [ "/bin/bash", "-c" ]

RUN apt-get -qq -y update && \
    DEBIAN_FRONTEND=noninteractive apt-get -qq -y install \
        R \
        python \
        libcairo2-dev \
        libtiff-dev \
        fftw3-dev \
        libharfbuzz-dev \
        libfribidi-dev && \
    apt-get -y autoclean && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt/lists/*

RUN git clone 

RUN Rscript -e 'install.packages("renv"); renv::'
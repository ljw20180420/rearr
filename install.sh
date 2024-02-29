#!/bin/bash

project_path="$(dirname $(realpath $0))"

# install dependencies
sudo apt-get -y update && \
sudo apt-get -y upgrade && \
sudo apt-get -y install libcairo2-dev libtiff-dev fftw3-dev libharfbuzz-dev libfribidi-dev imagemagick && \
sudo apt-get install pkg-config
sudo apt-get install build-essential gdb lcov pkg-config \
      libbz2-dev libffi-dev libgdbm-dev libgdbm-compat-dev liblzma-dev \
      libncurses5-dev libreadline6-dev libsqlite3-dev libssl-dev \
      lzma lzma-dev tk-dev uuid-dev zlib1g-dev
sudo apt-get -y autoremove

# install R
if [ ! -x "$project_path/R-4.3.2/bin/R" ]
then
    tar xzf "$project_path/R-4.3.2.tar.gz"
    cd "$project_path/R-4.3.2" || exit
    ./configure --enable-shared
    make
fi

# install python
if [ ! -x "$project_path/py312/bin/python3.12" ]
then
    tar xzf "$project_path/Python-3.12.1.tgz"
    cd "$project_path/Python-3.12.1" || exit
    ./configure --enable-optimizations --enable-shared --prefix="$project_path/py312" LDFLAGS="-Wl,-rpath $project_path/py312/lib"
    make
    make install
fi

# install pv
if [ ! -x "$project_path/pv-1.8.5/pv" ]
then
    tar xzf "$project_path/pv-1.8.5.tar.gz"
    cd "$project_path/pv-1.8.5" || exit
    ./configure
    make
fi

# install less
if [ ! -x "$project_path/less-643/less" ]
then
    unzip "$project_path/less-643.zip"
    cd "$project_path/less-643" || exit
    ./configure
    make
fi

# install kpLogo
if [ ! -x "$project_path/barcode/kpLogo-1.1/bin/kpLogo" ]
then
    tar xzf "$project_path/kpLogo-1.1.tar.gz" -C "$project_path/barcode"
    make -C "$project_path/barcode/kpLogo-1.1/src"
fi

cd $project_path

# install R packages
$project_path/R-4.3.2/bin/Rscript -e '
    if (!"pak" %in% rownames(installed.packages())) {
        install.packages("pak")
    }
    pak::pkg_install("languageserver")
    pak::pkg_install("jsonlite")
    pak::pkg_install("tidyverse")
    pak::pkg_install("ggforce")
    pak::pkg_install("waffle")
    pak::pkg_install("patchwork")
    pak::pkg_install("reticulate")
    pak::pkg_install("this.path")
    pak::pkg_install("ggtext")
    pak::pkg_install("Cairo")
    pak::pkg_install("omarwagih/ggseqlogo")
'

# install python packages
$project_path/py312/bin/pip3.12 install --upgrade pip
$project_path/py312/bin/pip3.12 install -r $project_path/requirements.txt

# simply append two paths to ~/.bashrc
if [[ $(grep -c "$(pwd)" ~/.bashrc) -eq 0 ]]
then
    printf "\n%s\n" "export PATH=\"$(pwd)\":"'$PATH' >> ~/.bashrc
fi

if [[ $(grep -c "$(pwd)"/barcode/tools ~/.bashrc) -eq 0 ]]
then
    printf "\n%s\n" "export PATH=\"$(pwd)/barcode/tools\":"'$PATH' >> ~/.bashrc
fi

# compile aligner
mkdir -p Rearrangement/build
cmake -DCMAKE_BUILD_TYPE=Release -S Rearrangement -B Rearrangement/build
make -C Rearrangement/build 
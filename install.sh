#!/bin/bash

project_path="$(dirname $(realpath $0))"

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

# install dependencies
sudo apt-get -y update && \
sudo apt-get -y upgrade && \
sudo apt-get -y install libcairo2-dev libtiff-dev fftw3-dev libharfbuzz-dev libfribidi-dev imagemagick && \
sudo apt-get -y autoremove

cd $project_path

# install R packages
Rscript --vanilla -e '
    args = commandArgs(trailingOnly = TRUE)
    renv::use_python(python = args[1])
    renv::restore()
' "$(which python3.10)"

# simply append three paths to ~/.bashrc
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
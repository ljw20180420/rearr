#!/bin/bash

#install dependencies
sudo apt-get -y update && \
sudo apt-get -y upgrade && \
sudo apt-get -y install libcairo2-dev libtiff-dev fftw3-dev libharfbuzz-dev libfribidi-dev

project_path="$(dirname $(realpath $0))"

cd "$project_path" || exit

# install R packages
Rscript --vanilla -e '
    args = commandArgs(trailingOnly = TRUE)
    renv::use_python(python = args[1])
    renv::restore()
' "$(which python3.10)"

# compile kpLogo for sgRNA analysis
make -C barcode/kpLogo/src

# simply append three paths to ~/.bashrc
if [[ $(grep -c "$(pwd)" ~/.bashrc) -eq 0 ]]
then
    printf "\n%s\n" "export PATH=\"$(pwd)\":"'$PATH' >> ~/.bashrc
fi

if [[ $(grep -c "$(pwd)"/barcode/tools ~/.bashrc) -eq 0 ]]
then
    printf "\n%s\n" "export PATH=\"$(pwd)/barcode/tools\":"'$PATH' >> ~/.bashrc
fi

mkdir -p Rearrangement/build
cmake -DCMAKE_BUILD_TYPE=Release -S Rearrangement -B Rearrangement/build
make -C Rearrangement/build

sudo ln -sf "$project_path/rearr_web.sh" /usr/local/bin
sudo ln -sf "$(which Rscript)" $project_path
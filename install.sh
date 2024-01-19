#!/bin/bash

cd "$(dirname "$0")" || exit

# install python packages
pip install -r requirements.txt

# install R packages
Rscript -e 'install.packages(c("tidyverse"), repos = "https://mirrors.sjtug.sjtu.edu.cn/cran/")'

# compile kpLogo for sgRNA analysis
make -C barcode/kpLogo/src

# simply append three paths to ~/.bashrc
if [[ $(grep -c "$(pwd)" ~/.bashrc) -eq 0 ]]
then
    printf "\n%s\n" "export PATH=\"$(pwd)\":"'$PATH' >> ~/.bashrc
fi

if [[ $(grep -c "$(pwd)"/tools ~/.bashrc) -eq 0 ]]
then
    printf "\n%s\n" "export PATH=\"$(pwd)/tools\":"'$PATH' >> ~/.bashrc
fi

if [[ $(grep -c "$(pwd)"/barcode/tools ~/.bashrc) -eq 0 ]]
then
    printf "\n%s\n" "export PATH=\"$(pwd)/barcode/tools\":"'$PATH' >> ~/.bashrc
fi

mkdir -p Rearrangement/build
cmake -DCMAKE_BUILD_TYPE=Release -S Rearrangement -B Rearrangement/build
make -C Rearrangement/build
if [[ $(grep -c "$(pwd)/Rearrangement/build" ~/.bashrc) -eq 0 ]]
then
    printf "\n%s\n" "export PATH=\"$(pwd)/Rearrangement/build\":"'$PATH' >> ~/.bashrc
fi

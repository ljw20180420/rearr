#!/bin/bash

cd $(dirname $0)

# install python packages
pip install -r requirements.txt

# simply append three paths to ~/.bashrc
if [[ $(grep -c $(pwd) ~/.bashrc) -eq 0 ]]
then
    printf "\n%s\n" "export PATH=\"$(pwd)\":"'$PATH' >> ~/.bashrc
fi
if [[ $(grep -c $(pwd)/tools ~/.bashrc) -eq 0 ]]
then
    printf "\n%s\n" "export PATH=\"$(pwd)/tools\":"'$PATH' >> ~/.bashrc
fi
if [[ $(grep -c $(pwd)/barcode/tools ~/.bashrc) -eq 0 ]]
then
    printf "\n%s\n" "export PATH=\"$(pwd)/barcode/tools\":"'$PATH' >> ~/.bashrc
fi
mkdir -p Rearrangement/build
cmake -DCMAKE_BUILD_TYPE=Release Rearrangement/build
make -C Rearrangement/build
if [[ $(grep -c $(pwd)/Rearrangement/build ~/.bashrc) -eq 0 ]]
then
    printf "\n%s\n" "export PATH=\"$(pwd)/Rearrangement/build\":"'$PATH' >> ~/.bashrc
fi
#!/bin/bash
# simply append three paths to ~/.bashrc

cd $(dirname $0)
if [[ $(grep -c $(pwd) ~/.bashrc) -eq 0 ]]
then
    printf "\n%s\n" "export PATH=\"$(pwd)\":"'$PATH' >> ~/.bashrc
fi
if [[ $(grep -c $(pwd)/tools ~/.bashrc) -eq 0 ]]
then
    printf "\n%s\n" "export PATH=\"$(pwd)/tools\":"'$PATH' >> ~/.bashrc
fi
cd Rearrangement
cmake -DCMAKE_BUILD_TYPE=Release .
make
if [[ $(grep -c $(pwd) ~/.bashrc) -eq 0 ]]
then
    printf "\n%s\n" "export PATH=\"$(pwd)\":"'$PATH' >> ~/.bashrc
fi
#!/bin/bash
cd $(dirname $0)
printf "\n%s\n" "export PATH=\"$(pwd)\":"'$PATH' >> ~/.bashrc
printf "\n%s\n" "export PATH=\"$(pwd)/tools\":"'$PATH' >> ~/.bashrc
cd Rearrangement
cmake -DCMAKE_BUILD_TYPE=Release .
make
printf "\n%s\n" "export PATH=\"$(pwd)\":"'$PATH' >> ~/.bashrc
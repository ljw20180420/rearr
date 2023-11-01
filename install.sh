#!/bin/bash
cd $(dirname $0)
printf "\n%s\n" "export PATH=\"$(pwd)\":"'$PATH' >> ~/.bashrc
printf "\n%s\n" "export PATH=\"$(pwd)/tools\":"'$PATH' >> ~/.bashrc
cd Rearrangement
cmake .
make
printf "\n%s\n" "export PATH=\"$(pwd)\":"'$PATH' >> ~/.bashrc
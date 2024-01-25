#!/bin/bash

cd "$(dirname "$0")" || exit

# install R packages
Rscript -e 'renv::restore()'


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

# bowtieindex="$HOME/hg19_with_bowtie2_index/hg19"
# read -p "path to bowtie2 index (default: /home/<user>/hg19_with_bowtie2_index/hg19):" inputbowtieindex
# test -n "$inputbowtieindex" && bowtieindex=$inputbowtieindex
# if [[ $(grep -c "export REARR_BOWTIE2_INDEX=" ~/.bashrc) -eq 0 ]]
# then
#     printf "\n%s\n" "export REARR_BOWTIE2_INDEX=\"$bowtieindex\"" >> ~/.bashrc
# fi
# genomeref="$HOME/hg19_with_bowtie2_index/hg19.fa"
# read -p "path to genome reference(default: /home/<user>/hg19_with_bowtie2_index/hg19.fa):" inputgenomeref
# test -n "$inputgenomeref" && genomeref=$inputbowtieindex
# if [[ $(grep -c "export REARR_GENOME_REF=" ~/.bashrc) -eq 0 ]]
# then
#     printf "\n%s\n" "export REARR_GENOME_REF=\"$genomeref\"" >> ~/.bashrc
# fi
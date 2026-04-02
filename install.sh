#!/bin/bash

# Change directory to the script directory.
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

default_prefix=${CONDA_PREFIX:-"${HOME}/.local"}
prefix=${prefix:-"${default_prefix}"}

# install rearrangement
core/Rearrangement/build.sh ${prefix}
mkdir -p ${prefix}/share/awk/
cp core/Rearrangement/correct_micro_homology.awk ${prefix}/share/awk/
# install removeDuplicates
cp core/removeDuplicates.sh ${prefix}/bin/
# install demultiplex
cp core/demultiplex/demultiplex.sh ${prefix}/bin/
cp core/demultiplex/getAlignPos.awk ${prefix}/share/awk/
# install runWorkFlow
cp ./runWorkFlow.sh ${prefix}/bin/
mkdir -p ${prefix}/share/rearr
cp ./workFlow.mak ${prefix}/share/rearr

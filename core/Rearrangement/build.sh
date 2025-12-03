#!/bin/bash

# Change directory to the script directory.
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

install_rearrangement=$1

rm -r build
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make

if [ "${install_rearrangement}" = "install" ]
then
    make install
fi

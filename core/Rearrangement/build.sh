#!/bin/bash

# Change directory to the script directory.
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

prefix=$1

rm -r build
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${prefix} ..
make
make install

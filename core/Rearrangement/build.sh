#!/bin/bash

# Change directory to the script directory.
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make

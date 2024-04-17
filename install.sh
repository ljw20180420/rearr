#!/bin/sh

# install rearrangement
mkdir -p Rearrangement/build
cd Rearrangement/build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
make install
cd -

# install pv
unzip pv-1.8.5
cd pv-1.8.5 
./configure LDFLAGS="-static"
make
make install
cd -

# install kpLogo
unzip kpLogo-1.1.zip
sed -i -r 's/(\$\(CC\) \$\(CFLAGS\))/\1 -static/; s/(gcc -O3)/\1 -static/' kpLogo-1.1/src/makefile
cd kpLogo-1.1/src
make
cp ../bin/kpLogo /usr/local/bin/
cd -

# install gawk
unzip gawk-5.3.0.zip
cd gawk-5.3.0
./configure LDFLAGS="-static"
make
make install
cd -

# install rearr_view.sh
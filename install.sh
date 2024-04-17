#!/bin/sh

for target in $@
do
    case $target in
        rearrangement|Rearrangement)
            # install rearrangement
            mkdir -p Rearrangement/build
            cd Rearrangement/build
            cmake -DCMAKE_BUILD_TYPE=Release ..
            make
            make install
            cd -
            ;;
        pv)
            # install pv
            unzip pv-1.8.5
            cd pv-1.8.5 
            ./configure LDFLAGS="-static"
            make
            make install
            cd -
            ;;
        kpLogo)
            # install kpLogo
            unzip kpLogo-1.1.zip
            sed -i -r 's/(\$\(CC\) \$\(CFLAGS\))/\1 -static/; s/(gcc -O3)/\1 -static/' kpLogo-1.1/src/makefile
            cd kpLogo-1.1/src
            make
            cp ../bin/kpLogo /usr/local/bin/
            cd -
            ;;
        gawk|awk)
            # install gawk
            unzip gawk-5.3.0.zip
            cd gawk-5.3.0
            ./configure LDFLAGS="-static"
            make
            make install
            cd -
            ;;
        correct)
            # install correct_micro_homology.AWK
            cp correct_micro_homology.AWK /usr/local/share/awk/
    esac
done
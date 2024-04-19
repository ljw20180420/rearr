#!/bin/sh

for target in $@
do
    case $target in
        rearrangement|Rearrangement)
            mkdir -p Rearrangement/build
            cd Rearrangement/build
            cmake -DCMAKE_BUILD_TYPE=Release ..
            make
            make install
            cd -
            ;;
        pv)
            unzip pv-1.8.5
            cd pv-1.8.5 
            ./configure LDFLAGS="-static"
            make
            make install
            cd -
            ;;
        kpLogo)
            unzip kpLogo-1.1.zip
            sed -i -r 's/(\$\(CC\) \$\(CFLAGS\))/\1 -static/; s/(gcc -O3)/\1 -static/' kpLogo-1.1/src/makefile
            cd kpLogo-1.1/src
            make
            cp ../bin/kpLogo /usr/local/bin/
            cd -
            ;;
        gawk|awk)
            unzip gawk-5.3.0.zip
            cd gawk-5.3.0
            ./configure LDFLAGS="-static"
            make
            make install
            cd -
            ;;
        correct)
            cp correct_micro_homology.awk /usr/local/share/awk/
            ;;
        demultiplex)
            cp demultiplex/demultiplex.sh /usr/local/bin/
            cp demultiplex/endOfSpliterPos.awk /usr/local/share/awk/
            cp demultiplex/cumulateToMapCutAdaptSpliter.awk /usr/local/share/awk/
            ;;
        sx)
            # install getSxCsvFileRef
            cp sx/getSxCsvFileRef/getSxCsvFileRef.sh /usr/local/bin/
            cp sx/getSxCsvFileRef/getSxCsvFileTarget.pl /usr/local/bin/
            cp sx/getSxCsvFileRef/getSxRefFile.pl /usr/local/bin/
            cp sx/getSxCsvFileRef/sxTargetSam2Bed.awk /usr/local/share/awk/
            # install sxInderSpliter
            cp sx/sxIndexSpliter.sh /usr/local/bin/
            # install sxInferCsvfileFromFastq
            cp sx/sxInferCsvfileFromFastq.sh /usr/local/bin/
            ;;
    esac
done
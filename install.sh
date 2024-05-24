#!/bin/bash
# Usage: sudo ./install.sh core sx

for target in $@
do
    case $target in
        apt|core)
            apt-get update && apt-get install -y --no-install-recommends unzip build-essential libncurses5-dev gawk bowtie2 cutadapt samtools cmake bedtools
            ;;&
        rearrangement|Rearrangement|core)
            mkdir -p Rearrangement/build
            cd Rearrangement/build
            cmake -DCMAKE_BUILD_TYPE=Release ..
            make
            make install
            cd -
            ;;&
        pv|core)
            unzip -o pv-1.8.5
            cd pv-1.8.5 
            ./configure LDFLAGS="-static"
            make
            make install
            cd -
            ;;&
        kpLogo)
            unzip -o kpLogo-1.1.zip
            sed -i -r 's/(\$\(CC\) \$\(CFLAGS\))/\1 -static/; s/(gcc -O3)/\1 -static/' kpLogo-1.1/src/makefile
            cd kpLogo-1.1/src
            make
            cp ../bin/kpLogo /usr/local/bin/
            cd -
            ;;&
        correct|core)
            cp correct_micro_homology.awk /usr/share/awk/
            ;;&
        removeDup|core)
            cp removeDuplicates.sh /usr/local/bin/
            ;;&
        demultiplex|core)
            cp demultiplex/demultiplex.sh /usr/local/bin/
            cp demultiplex/endOfSpliterPos.awk /usr/share/awk/
            ;;&
        sx)
            # install getSxCsvFileRef
            cp sx/getSxCsvFileRef/getSxCsvFileRef.sh /usr/local/bin/
            cp sx/getSxCsvFileRef/getSxCsvFileTarget.pl /usr/local/bin/
            cp sx/getSxCsvFileRef/getSxRefFile.pl /usr/local/bin/
            cp sx/getSxCsvFileRef/sxTargetSam2Bed.awk /usr/share/awk/
            # install cutR2Adapter
            cp sx/sxCutR2AdapterFilterCumulate/sxCutR2AdapterFilterCumulate.sh /usr/local/bin/
            cp sx/sxCutR2AdapterFilterCumulate/sxCumulateToMapCutAdaptSpliter.awk /usr/share/awk/
            # install sxInderSpliter
            cp sx/sxExtractSpliter.sh /usr/local/bin/
            ;;&
    esac
done

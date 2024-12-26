#!/bin/bash
# Usage: sudo ./install.md core sx

for target in $@
do
    case $target in
        apt|core)
            apt-get update && apt-get install -y --no-install-recommends unzip build-essential libncurses5-dev gawk bowtie2 cutadapt samtools cmake bedtools
            ;;&
        rearrangement|Rearrangement|core)
            mkdir -p core/Rearrangement/build
            cd core/Rearrangement/build
            cmake -DCMAKE_BUILD_TYPE=Release ..
            make
            make install
            cd -
            cp core/Rearrangement/correct_micro_homology.awk /usr/share/awk/
            ;;&
        removeDup|core)
            cp core/removeDuplicates.md /usr/local/bin/
            ;;&
        demultiplex|core)
            cp core/demultiplex/demultiplex.md /usr/local/bin/
            cp core/demultiplex/getAlignPos.awk /usr/share/awk/
            ;;&
        sx)
            # install getSxCsvFileRef
            cp sx/getSxCsvFileRef/getSxCsvFileRef.md /usr/local/bin/
            cp sx/getSxCsvFileRef/getSxCsvFileTarget.pl /usr/local/bin/
            cp sx/getSxCsvFileRef/getSxRefFile.pl /usr/local/bin/
            cp sx/getSxCsvFileRef/sxTargetSam2Bed.awk /usr/share/awk/
            # install cutR2Adapter
            cp sx/sxCutR2AdapterFilterCumulate/sxCutR2AdapterFilterCumulate.md /usr/local/bin/
            cp sx/sxCutR2AdapterFilterCumulate/sxCumulateToMapCutAdaptSpliter.awk /usr/share/awk/
            # install sxInderSpliter
            cp sx/sxExtractSpliter.md /usr/local/bin/
            ;;&
    esac
done

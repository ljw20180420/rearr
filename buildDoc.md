#!/bin/bash

add_header()
{
    title=$1
    permalink=$2
    sed "1i ---\ntitle: \"$title\"\npermalink: $permalink\ntoc: true\n---\n"
}

wrap_script()
{
    script_type=$1
    sed -e "1i ~~~$script_type" -e '$a~~~'
}

add_header "Remove duplicates" "/core/remove-duplicates/" < core/removeDuplicates.md > docs/_docs/removeDuplicates.md

add_header "Demultiplex" "/core/demultiplex/" < core/demultiplex/demultiplex.md > docs/_docs/demultiplex.md
wrap_script awk < core/demultiplex/getAlignPos.awk | add_header "Get alignment position" "/core/demultiplex/get-alignment-position/" > docs/_docs/getAlignPos.md

add_header "Chimeric alignment" "/core/rearr/" < core/Rearrangement/rearr.md > docs/_docs/rearr.md
wrap_script awk < core/Rearrangement/correct_micro_homology.awk | add_header "Correct micro-homology" "/core/rearr/correct-micro-homology/" > docs/_docs/correct_micro_homology.md
cd core/Rearrangement
doxygen
cd -

add_header "Shi Xing extract spliter" "/sx/sx-extract-spliter/" < sx/sxExtractSpliter.md > docs/_docs/sxExtractSpliter.md

add_header "Shi Xing post-process from demultiplex to rearr" "/sx/sx-cut-r2-adapter-filter-cumulate/" < sx/sxCutR2AdapterFilterCumulate/sxCutR2AdapterFilterCumulate.md > docs/_docs/sxCutR2AdapterFilterCumulate.md
wrap_script awk < sx/sxCutR2AdapterFilterCumulate/sxCumulateToMapCutAdaptSpliter.awk | add_header "Accumulate adjacent duplicated queries" "/sx/sx-cut-r2-adapter-filter-cumulate/sx-cumulate-to-map-cut-adapt-spliter/" > docs/_docs/sxCumulateToMapCutAdaptSpliter.md

add_header "get Shi Xing csvfile reference" "/sx/get-sx-csvfile-ref/" < sx/getSxCsvFileRef/getSxCsvFileRef.md > docs/_docs/getSxCsvFileRef.md
wrap_script perl < sx/getSxCsvFileRef/getSxCsvFileTarget.pl | add_header "get Shi Xing csvfile target" "/sx/get-sx-csvfile-ref/get-sx-csvfile-target/" > docs/_docs/getSxCsvFileTarget.md
wrap_script awk < sx/getSxCsvFileRef/sxTargetSam2Bed.awk | add_header "Shi Xing target sam to bed" "/sx/get-sx-csvfile-ref/sx-target-sam-2-bed/" > docs/_docs/sxTargetSam2Bed.md
wrap_script perl < sx/getSxCsvFileRef/getSxRefFile.pl | add_header "get Shi Xing reference file" "/sx/get-sx-csvfile-ref/get-sx-ref-file/" > docs/_docs/getSxRefFile.md

cd docs
bundle exec jekyll serve

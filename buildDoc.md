#!/bin/bash

sed '1i ---\ntitle: "Remove duplicates"\npermalink: /remove-duplicates/\ntoc: true\n---\n' core/removeDuplicates.md > docs/_docs/removeDuplicates.md

sed '1i ---\ntitle: "Demultiplex"\npermalink: /demultiplex/\ntoc: true\n---\n' core/demultiplex/demultiplex.md > docs/_docs/demultiplex.md
sed -e '1i ---\ntitle: "Get alignment position"\npermalink: /demultiplex/get-alignment-position/\ntoc: true\n---\n\n~~~awk' -e '$a~~~' core/demultiplex/getAlignPos.awk > docs/_docs/getAlignPos.md

cp core/Rearrangement/rearr.md docs/_docs/
sed -e '1i ---\ntitle: "Correct micro-homology"\npermalink: /rearr/correct-micro-homology/\ntoc: true\n---\n\n~~~awk' -e '$a~~~' core/Rearrangement/correct_micro_homology.awk > docs/_docs/correct_micro_homology.md
cd core/Rearrangement
doxygen
cd -

cd docs
bundle exec jekyll serve
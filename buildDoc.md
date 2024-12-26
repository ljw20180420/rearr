#!/bin/bash

sed '1i ---\ntitle: "Remove duplicates"\npermalink: /remove-duplicates/\ntoc: true\n---\n' pre-post-process/removeDuplicates.md > docs/_docs/removeDuplicates.md

cd Rearrangement
doxygen
cp rearr.md ../docs/_docs/

cd ../docs
bundle exec jekyll serve
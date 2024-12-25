#!/bin/bash

cd Rearrangement
doxygen
cp rearr.md ../docs/_docs/
cd ../docs
bundle exec jekyll serve
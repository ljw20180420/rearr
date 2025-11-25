#!/bin/bash

# change to the dir of the script
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

manim render --save_sections main.py MainScene
manim render --save_sections algorithm.py Forward

#!/bin/bash

# change to the dir of the script
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

manim render --save_sections main.py MainScene
manim render --save_sections algorithm.py Forward

./render.sh
cp \
    manim/media/videos/main/480p30/MainScene.mp4 \
    manim/media/videos/algorithm/480p30/Forward.mp4 \
    docs/assets/videos/

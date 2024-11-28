#!/bin/bash

# This is a script to run before docker compose

chmod a+w docker-images/shiny-server/diffloopAnalysis/www
chmod a+w docker-images/shiny-server/diffloopAnalysisPair/www
chmod a+w docker-images/shiny-server/downstreamAnalysis/www
mkdir -p flask_project/tmp
chmod a+w flask_project/tmp
mkdir -p chat-ui/data/db
chmod a+w chat-ui/data/db
mkdir -p chat-ui/llama.cpp
wget https://hf-mirror.com/microsoft/Phi-3-mini-4k-instruct-gguf/resolve/main/Phi-3-mini-4k-instruct-q4.gguf -O chat-ui/llama.cpp/Phi-3-mini-4k-instruct-q4.gguf

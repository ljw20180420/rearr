#!/bin/bash

# This is a script to run before docker compose

chmod a+w docker-images/shiny/logs
chmod a+w docker-images/shiny/apps/diffloopAnalysis/www
chmod a+w docker-images/shiny/apps/diffloopAnalysisPair/www
chmod a+w docker-images/shiny/apps/downstreamAnalysis/www
chmod a+w docker-images/flask/flask_project/tmp
chmod a+w docker-images/chat-ui/data/db
if ! [ -e "docker-images/chat-ui/llama.cpp/Phi-3-mini-4k-instruct-q4.gguf" ]
then
    wget https://hf-mirror.com/microsoft/Phi-3-mini-4k-instruct-gguf/resolve/main/Phi-3-mini-4k-instruct-q4.gguf -O docker-images/chat-ui/llama.cpp/Phi-3-mini-4k-instruct-q4.gguf
fi
docker compose down
docker compose up -d
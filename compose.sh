#!/bin/bash

# This is a script to run before docker compose

chmod a+w docker-images/shiny/logs
chmod a+w docker-images/shiny/apps/diffloopAnalysis/www
chmod a+w docker-images/shiny/apps/diffloopAnalysisPair/www
chmod a+w docker-images/shiny/apps/downstreamAnalysis/www
chmod a+w docker-images/flask/flask_project/tmp
chmod a+w docker-images/chat-ui/chat
chmod a+w docker-images/chat-ui/TEI
chmod a+w docker-images/chat-ui/TGI

docker compose down

# docker compose build

# docker compose push

# docker compose pull

# for file in $(ls docker-images/chat-ui/data/db); do rm -rf "docker-images/chat-ui/data/db/$file"; done

docker compose up -d

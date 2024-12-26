#!/bin/bash

# This is a script to run before docker compose

# shiny
chmod a+w docker-images/shiny/logs
chmod a+w docker-images/shiny/apps/diffloopAnalysis/www
chmod a+w docker-images/shiny/apps/diffloopAnalysisPair/www
chmod a+w docker-images/shiny/apps/downstreamAnalysis/www

# workflow
chmod a+w docker-images/flask/flask_project/tmp

docker compose down

docker compose build

docker compose push

docker compose pull

docker compose up -d

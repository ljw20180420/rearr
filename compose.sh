#!/bin/bash

# change to the dir of the script
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# This is the start script for docker compose.

# build vue3
cd docker-images/flask/vue_project
npm run build
cd - 

# shiny
chmod a+w docker-images/shiny/logs
chmod a+w docker-images/shiny/apps/downstreamAnalysis/www

# workflow
chmod a+w docker-images/flask/tmp

# Stop and remove previous up
docker compose down
# Pull remote images
docker compose pull
# Up in the background
docker compose up -d

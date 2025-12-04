#!/bin/bash

# change to the dir of the script
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# This is the start script for docker compose.

# build vue3
cd flask_app/vue_project
npm run build
cd - 

# shiny
chmod a+w shiny/logs
chmod a+w shiny/apps/downstreamAnalysis/www

# workflow
chmod a+w flask_app/tmp

# Stop and remove previous up
docker compose down
# Pull remote images
docker compose pull
# Up in the background
docker compose up -d

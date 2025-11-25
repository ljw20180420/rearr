#!/bin/bash

# change to the dir of the script
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# This is the start script for docker compose. It handles volume permissions, stop and remove previous up if the compose is already running, pull the latest remote images, and finally up the compose in the background.

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

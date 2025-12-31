#!/bin/bash

# change to the dir of the script
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# This is the start script for docker compose.

# build jekyll
pushd jekyll
bundle exec jekyll build
popd

# build vue3
pushd flask_app/vue_project
npm run build
popd

# shiny
mkdir -p shiny/docker_apps/downstreamAnalysis
cp shiny/apps/downstreamAnalysis/app.R shiny/docker_apps/downstreamAnalysis/app.R
cp -r shiny/apps/downstreamAnalysis/library shiny/docker_apps/downstreamAnalysis/
chmod a+w shiny/docker_apps/downstreamAnalysis

# Stop and remove previous up
docker compose down
# Pull remote images
docker compose pull
# Up in the background
docker compose up -d

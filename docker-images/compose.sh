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
mkdir -p docker_mounts/shiny_apps/downstreamAnalysis
cp shiny/apps/downstreamAnalysis/app.R docker_mounts/shiny_apps/downstreamAnalysis/app.R
cp -r shiny/apps/downstreamAnalysis/library docker_mounts/shiny_apps/downstreamAnalysis/
chmod a+w docker_mounts/shiny_apps/downstreamAnalysis

# flask
mkdir -p docker_mounts/flask_app/vue_project
cp flask_app/__init__.py flask_app/__main__.py flask_app/tasks.py docker_mounts/flask_app/
cp -r flask_app/vue_project/dist docker_mounts/flask_app/vue_project/
chmod a+w docker_mounts/flask_app

# Stop and remove previous up
docker compose down
# Pull remote images
docker compose pull
# Up in the background
docker compose up -d

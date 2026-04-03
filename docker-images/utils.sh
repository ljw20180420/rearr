#!/bin/bash

build_jekyll() {
    pushd jekyll
    bundle update --bundler
    bundle lock --update
    bundle install
    bundle exec jekyll build
    popd
}

build_shiny() {
    mkdir -p docker_mounts/shiny_apps/downstreamAnalysis
    cp shiny/apps/downstreamAnalysis/app.R docker_mounts/shiny_apps/downstreamAnalysis/app.R
    cp -r shiny/apps/downstreamAnalysis/library docker_mounts/shiny_apps/downstreamAnalysis/
    chmod a+w docker_mounts/shiny_apps/downstreamAnalysis
}

build_vue() {
    pushd flask_app/vue_project
    npm install
    npm update
    npm run build
    popd
}

build_flask() {
    build_vue
    mkdir -p docker_mounts/flask_app/vue_project
    cp flask_app/__init__.py flask_app/__main__.py flask_app/tasks.py docker_mounts/flask_app/
    cp -r flask_app/vue_project/dist docker_mounts/flask_app/vue_project/
    chmod a+w docker_mounts/flask_app
}

build_webUI() {
    build_jekyll
    build_shiny
    build_flask
}
#!/bin/bash

# This is a script to run before docker compose

mkdir -p flask_project/tmp
chmod a+w flask_project/tmp
mkdir -p shinyApps/downstreamAnalysis/www
chmod a+w shinyApps/downstreamAnalysis/www
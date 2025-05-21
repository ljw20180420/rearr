#!/bin/bash

shopt -s expand_aliases

alias ~~~=":<<'~~~bash'"

:<<'~~~bash'

# Usage
```bash
./compose.md
```

# Introduction
This is the start script for docker compose. It handles volume permissions, stop and remove previous up if the compose is already running, pull the latest remote images, and finally up the compose in the background.

# Source
~~~bash
# shiny
chmod a+w docker-images/shiny/logs
chmod a+w docker-images/shiny/apps/diffloopAnalysis/www
chmod a+w docker-images/shiny/apps/diffloopAnalysisPair/www
chmod a+w docker-images/shiny/apps/downstreamAnalysis/www

# workflow
chmod a+w docker-images/flask/tmp

# Stop and remove previous up
docker compose down
# Pull remote images
docker compose pull
# Up in the background
docker compose up -d
~~~

~~~bash
alias ~~~=":" # This suppresses a warning and is not part of source.
~~~

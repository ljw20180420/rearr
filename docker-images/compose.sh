#!/bin/bash

# change to the dir of the script
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# This is the start script for docker compose.

. utils.sh
build_webUI

# Stop and remove previous up
docker compose down
# Pull remote images
docker compose pull
# Up in the background
docker compose up -d

#!/bin/bash

# change to the dir of the script
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

. utils.sh
build_flask

# change to the dir of docker_mounts
cd docker_mounts

docker run --rm -d -p 5672:5672 rabbitmq:latest
docker run --rm -d -p 6379:6379 redis:latest

celery -A flask_app.tasks \
    -b amqp://localhost:5672 \
    --result-backend redis://localhost:6379/0 \
    worker --loglevel=INFO --detach --concurrency 1

# port is 8080
python -m flask_app \
    -b amqp://localhost:5672 \
    --result-backend redis://localhost:6379/0

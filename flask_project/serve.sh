#!/bin/bash

systemctl restart rabbitmq-server.service
systemctl restart celery.service

# flaskPath=$(dirname $(realpath $0))
# cd flaskPath

# cd ${flaskPath}/vue_project
# npm run build
# cd ${flaskPath}
# celery -A rearr_web worker --loglevel="INFO" --concurrency=2 &
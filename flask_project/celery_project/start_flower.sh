#!/bin/bash

flaskPath=$(dirname $(realpath $0))/..
docker run -d --restart=always -v ${flaskPath}:/data -p 5555:5555 mher/flower celery --app=celery_project.app flower
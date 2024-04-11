#!/bin/bash

if ! docker ps | grep "celery_flower"
then
    if docker ps -a | grep "celery_flower"
    then
        docker container start celery_flower
    else
        flaskPath=$(dirname $(realpath $0))/..
        docker run --name celery_flower -d --restart=always -v ${flaskPath}:/data -p 5555:5555 mher/flower celery --app=celery_project.app flower
    fi
fi



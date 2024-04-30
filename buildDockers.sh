#!/bin/bash
# It is suggested that not build images in compose.yaml, so I build it here.
docker build --progress=plain -f celery_worker.df -t my_celery_worker .
docker build --progress=plain -f celery_flower.df -t my_celery_flower .
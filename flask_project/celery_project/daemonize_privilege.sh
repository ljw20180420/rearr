#!/bin/bash

flaskPath=$(dirname $(dirname $(realpath $0)))
cp ${flaskPath}/celery_project/celery.service /etc/systemd/system/
systemctl daemon-reload
systemctl restart celery.service
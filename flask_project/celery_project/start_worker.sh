#!/bin/bash

projectPath=$(dirname $(realpath $0))/../..
# Generally, there is no point in running more than one worker on a particular machine unless you want to do routing.
cd "${projectPath}/flask_project"
celery -A celery_project.app worker --loglevel="INFO" --concurrency=2
#!/bin/bash

flaskPath=$(dirname $(dirname $(realpath $0)))
User=$(whoami)
Group=$(id -gn)
WorkingDirectory=${flaskPath}
ExecStart="$(which celery) -A celery_project.app multi start w1 --pidfile=${WorkingDirectory}/run/%n.pid --logfile=${WorkingDirectory}/log/%n%I.log --loglevel=INFO"
ExecStop="$(which celery) multi stopwait w1 --pidfile=${WorkingDirectory}/run/%n.pid --logfile=${WorkingDirectory}/log/%n%I.log --loglevel=INFO"
ExecReload=$(sed -r 's/\sstart\s/restart/' <<<${ExecStart})
sed -e 's|^User=|User='${User}'|;' \
    -e 's|^Group=|Group='${Group}'|;' \
    -e 's|^WorkingDirectory=|WorkingDirectory='${WorkingDirectory}'|;' \
    -e 's|^ExecStart=|ExecStart='"${ExecStart}"'|;' \
    -e 's|^ExecStop=|ExecStop='"${ExecStop}"'|;' \
    -e 's|^ExecReload=|ExecReload='"${ExecReload}"'|;' \
    ${flaskPath}/celery_project/celery.service.template >${flaskPath}/celery_project/celery.service
    
    

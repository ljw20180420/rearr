#!/bin/bash
# Usage
# rearr_render_docker.sh input.alg.cut1.cut2.ref1len

algfile=$1
if [ ! -s $algfile ]
then
    echo "$algfile does not exist or is empty"
    exit 1
fi

project_path="$(dirname $(realpath $0))"
inputpath="$(dirname $(realpath $algfile))"
docker run --rm --mount type=bind,src="${inputpath}",dst="/app/data" rearr:auto ./rearr_render.sh ./data/$(basename $algfile)
# docker run -it --rm --mount type=bind,src="${inputpath}",dst="/app/data" rearr:auto /bin/bash
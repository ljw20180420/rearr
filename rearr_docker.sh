#!/bin/bash
# usage: rearr_docker.sh relative_path_to_shell_script args
# relative means relative to rearr_docker.sh (not the current working directory)
# the first arg in args is an input file, the left args are parameters
shell_script=$1
inputfile=$2
shift 2
inputpath=$(dirname $(realpath $inputfile))

if [[ $shell_script == *rearr_view.sh ]]
then
    docker run -it --rm --mount type=bind,src="$inputpath",dst="/app/data" rearr:auto ./$shell_script ./data/$(basename $inputfile) $@
else
    docker run --rm --mount type=bind,src="$inputpath",dst="/app/data" rearr:auto ./$shell_script ./data/$(basename $inputfile) $@
fi
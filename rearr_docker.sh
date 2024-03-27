#!/bin/bash

project_path=$(dirname $(realpath $0))
exec=./$(realpath --relative-base=${project_path} $(find ${project_path} -name "$1"))
shift
files=$(echo $@ | gawk -F "--" '{print $1}')
vars=$(echo $@ | gawk -F "--" '{print $2}')
declare -a srcs
for file in $files
do
    file=$(realpath ${file})
    filedir=$(dirname ${file})
    insrcs=0
    for src in "${srcs[@]}"
    do
        if [ $src == ${filedir} ]
        then
            insrcs=1
            break
        fi
    done
    if [ $insrcs == 0 ]
    then
        srcs[${#srcs[@]}]=${filedir}
    fi
    fileDocker="${fileDocker} ${file}"
done
for src in "${srcs[@]}"
do
    mount="${mount} --mount type=bind,src=${src},dst=${src}"
done
docker run --rm $mount rearr:auto $exec $fileDocker $vars


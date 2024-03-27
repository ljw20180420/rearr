#!/bin/bash
project_path="$(dirname $(realpath $0))"
# build docker
docker bulid --progress plain -t rearr:auto .
# index csvfiles
read -ep "path to genome reference:" genomeref
"${project_path}/rearr_docker.sh" index_spliter.sh $genomeref $(ls ${project_path}/barcode/sx/csvfiles/*.csv)
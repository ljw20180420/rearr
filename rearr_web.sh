#!/bin/bash
fqfile=$PATH_INFO
read ref1 ref2 cut1 cut2 NGGCCNtype1 NGGCCNtype2 < <(sed 's/&/ /g' <<<$QUERY_STRING)
printf "%s<br>" "fqfile=$fqfile" "ref1=$ref1" "ref2=$ref2" "cut1=$cut1" "cut2=$cut2" "NGGCCNtype1=$NGGCCNtype1" "NGGCCNtype2=$NGGCCNtype2"
fqfile="/var/www/rearrangement.org/public_html/$fqfile"
project_path="$(dirname $(realpath $0))"
$project_path/rearr_run.sh $fqfile $ref1 $ref2 $cut1 $cut2 $NGGCCNtype1 $NGGCCNtype2
export HOME=~"$(cut -d/ -f3 $project_path)"
echo $HOME
$project_path/rearr_render.sh $fqfile.alg.$cut1.$cut2.${#ref1}
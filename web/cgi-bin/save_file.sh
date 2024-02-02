#!/bin/bash
printf "Content-Type: text\n\n"

project_path="$(dirname $(realpath $0))/../.."
cat >"$project_path/web/jobs/tmp.fq"
#!/bin/bash
project_path="$(dirname $(realpath $0))"
$project_path/py312/bin/python3.12 -m http.server --cgi --directory $project_path/web --bind localhost 8000
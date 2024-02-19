#!/bin/bash
project_path="$(dirname $(realpath $0))"
$project_path/py312/bin/python3.12 -m flask --app $project_path/web/rearr_web.py run -h localhost -p 8000
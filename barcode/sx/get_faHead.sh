#!/bin/bash
# Usage: get_faHead.sh <csvfile
cut -d, -f1 | sed 's/^/>BC_/'
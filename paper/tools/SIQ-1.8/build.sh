#!/bin/bash

# change to the dir of the script
cd $(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

rm -rf SIQ/bin
# -g Generate all debugging info
javac \
	-classpath 'SIQ/libraries/*' \
	-sourcepath SIQ/src \
	-d SIQ/bin \
	SIQ/src/my/my_SIQ.java

#!/bin/bash

shopt -s expand_aliases

alias ~~~=":<<'~~~bash'"

:<<'~~~bash'

# Usage
```bash
./loginWorker.md
```

# Introduction
Login to a temporary container of the image `ghcr.io/ljw20180420/celery_worker`. The container will be remove automatically after exiting. This script accept a `dataPath` parameter (default to `./data`). The `dataPath` contain the input `.(fastq|fq)[.(gz|zip)]` files. It maps to `/app/data` as a volume. The `dataPath` may optionally contain the genome and the corresponding bowtie2 index. These are necessary for using the automatic reference extraction of the in-house sx module. The `WorkingDir` of `ghcr.io/ljw20180420/celery_worker` is `/app`.

# Source
~~~bash
read -ep "please input the full path of directory containing .(fq|fastq)[.(zip|gz)] inputs and other necessary data:" dataPath
if [ -z $dataPath ]
then
    echo "No dataPath is specified. Use default ./data as dataPath." >&2
    dataPath='./data'
fi
docker run -it --rm \
-v "./test:/app/test" \
-v "$dataPath:/app/data" \
ghcr.io/ljw20180420/celery_worker /bin/bash
~~~

~~~bash
alias ~~~=":" # This suppresses a warning and is not part of source.
~~~

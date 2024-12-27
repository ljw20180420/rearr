#!/bin/bash

shopt -s expand_aliases

alias ~~~=":<<'~~~bash'"

:<<'~~~bash'

# Usage
```bash
./loginWorker.md
```

# Introduction
Login to a temporary container of the image `ghcr.io/ljw20180420/celery_worker`. The container will be remove automatically after exiting. Both a genome and the corresponding bowtie2 index is necessary for using the automatic reference extraction of the in-house sx module. You can skip (press `enter`) both if you do not use automatic reference extraction. This script also require a mandatory data path. The data path contain the input `(fastq|fq).(gz|zip)` files. It maps to `/app/data` as a volume. The `WorkingDir` of `ghcr.io/ljw20180420/celery_worker` is `/app`.

# Source
~~~bash
# select genome path
if [ ! -e genome/genome.fa ] || [ ! -e genome/genome.1.bt2 ] || [ ! -e genome/genome.2.bt2 ] || [ ! -e genome/genome.3.bt2 ] || [ ! -e genome/genome.4.bt2 ] || [ ! -e genome/genome.rev.1.bt2 ] || [ ! -e genome/genome.rev.2.bt2 ]
then
    mkdir -p genome
    # example: /home/abc/hg19/hg19.fa
    read -ep "please input the full path of genome fasta file:" genome
    if [ -z $genome ]
    then
        echo "No genome is provided. Reference extraction of the in-house sx module will not be usable."
    else
        ln $genome ./genome/genome.fa
    fi
    # example: /home/abc/hg19/hg19
    read -ep "please input the full path of genome bowtie2 index:" bowtie2index
    if [ -z $bowtie2index ]
    then
        echo "No genome index is provided. Reference extraction of the in-house sx module will not be usable."
    else
        for ext in 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2
        do
            ln $bowtie2index.$ext ./genome/genome.$ext
        done
    fi
    chown $(id -un):$(id -gn) genome
fi

read -ep "please input the full path of directory containing .fastq\.fq(.gz):" dataPath
while [ -z $dataPath ]
do
    echo "must specify dataPath" >&2
    read -ep "please input the full path of directory containing .fastq\.fq(.gz):" dataPath
done
docker run -it --rm \
-v "./genome:/app/genome" \
-v "./sx/csvfiles:/app/sx/csvfiles" \
-v "./test:/app/test" \
-v "$dataPath:/app/data" \
ghcr.io/ljw20180420/celery_worker /bin/bash
~~~

~~~bash
alias ~~~=":" # This suppresses a warning and is not part of source.
~~~

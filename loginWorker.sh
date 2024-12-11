#!/bin/bash
# Usage: loginWorker.sh

# select genome path
if [ ! -e genome/genome.fa ] || [ ! -e genome/genome.1.bt2 ] || [ ! -e genome/genome.2.bt2 ] || [ ! -e genome/genome.3.bt2 ] || [ ! -e genome/genome.4.bt2 ] || [ ! -e genome/genome.rev.1.bt2 ] || [ ! -e genome/genome.rev.2.bt2 ]
then
    mkdir -p genome
    # example: /home/abc/hg19/hg19.fa
    read -ep "please input the full path of genome fasta file:" genome
    ln $genome ./genome/genome.fa
    # example: /home/abc/hg19/hg19
    read -ep "please input the full path of genome bowtie2 index:" bowtie2index
    for ext in 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2
    do
        ln $bowtie2index.$ext ./genome/genome.$ext
    done
    chown $(id -un):$(id -gn) genome
fi

read -ep "please input the full path of directory containing .fastq\.fq(.gz):" dataPath
while [ -z $dataPath ]
do
    echo "must specify dataPath" >&2
    read -ep "please input the full path of directory containing .fastq\.fq(.gz):" dataPath
fi
docker run -it --rm \
-v "./genome:/app/genome" \
-v "./sx/csvfiles:/app/sx/csvfiles" \
-v "./test:/app/test" \
-v "$dataPath:/app/data" \
ljwdocker1989/celery_worker /bin/bash

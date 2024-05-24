#!/bin/bash

docker run -it --rm -v "./hg19:/app/hg19" -v "./sx/csvfiles:/app/sx/csvfiles" -v "./test:/app/test" ljwdocker1989/celery_worker /bin/bash
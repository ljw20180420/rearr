#!/bin/bash
ref1=$(sed -n '2p' /home/ljw/wuqiang/sx/sx_lcy/test/ref12.fa)
ref2=$(sed -n '4p' /home/ljw/wuqiang/sx/sx_lcy/test/ref12.fa)
rearr_run.sh /home/ljw/wuqiang/sx/sx_lcy/test/random.fq $ref1 $ref2 100 100 NGG NGG
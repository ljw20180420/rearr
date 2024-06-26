#!/bin/bash

# docker run -it --rm -v ./certbot/letsencrypt/:/etc/letsencrypt/ -v ./certbot/lib/letsencrypt/:/var/lib/letsencrypt/ certbot/certbot certonly -d qiangwulab.sjtu.edu.cn -d www.qiangwulab.sjtu.edu.cn -d gallery.qiangwulab.sjtu.edu.cn -d workflow.qiangwulab.sjtu.edu.cn -d flower.qiangwulab.sjtu.edu.cn -d shiny.qiangwulab.sjtu.edu.cn -d rearr.xyz -d www.rearr.xyz -d gallery.rearr.xyz -d workflow.rearr.xyz -d flower.rearr.xyz -d shiny.rearr.xyz --standalone --preferred-challenges http --agree-tos -m ljw2017@sjtu.edu.cn --cert-name wulab

# docker run -it --rm -v ./certbot/letsencrypt/:/etc/letsencrypt/ -v ./certbot/lib/letsencrypt/:/var/lib/letsencrypt/ certbot/certbot certonly -d qiangwulab.sjtu.edu.cn -d *.qiangwulab.sjtu.edu.cn -d rearr.xyz -d *.rearr.xyz --standalone --preferred-challenges http --agree-tos -m ljw2017@sjtu.edu.cn --cert-name wulab

docker run -it --rm -v ./certbot/letsencrypt/:/etc/letsencrypt/ -v ./certbot/lib/letsencrypt/:/var/lib/letsencrypt/ certbot/certbot certonly -d rearr.xyz -d *.rearr.xyz --standalone --preferred-challenges http --agree-tos -m ljw2017@sjtu.edu.cn --cert-name wulab

# if ! grep -E "0 12 \* \* \* docker run -it --rm -v ./certbot/letsencrypt/:/etc/letsencrypt/ -v ./certbot/lib/letsencrypt/:/var/lib/letsencrypt/ certbot/certbot renew -q" /etc/crontab
# then
#     echo "0 12 * * * docker run -it --rm -v ./certbot/letsencrypt/:/etc/letsencrypt/ -v ./certbot/lib/letsencrypt/:/var/lib/letsencrypt/ certbot/certbot renew -q" >>/etc/crontab
# fi
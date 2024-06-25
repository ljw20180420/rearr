#!/bin/bash

docker run -it --rm -v ./certbot/letsencrypt/:/etc/letsencrypt/ -v ./certbot/lib/letsencrypt/:/var/lib/letsencrypt/ certbot/certbot certonly -d qiangwulab.sjtu.edu.cn -d www.qiangwulab.sjtu.edu.cn -d gallery.qiangwulab.sjtu.edu.cn -d workflow.qiangwulab.sjtu.edu.cn -d flower.qiangwulab.sjtu.edu.cn -d shiny.qiangwulab.sjtu.edu.cn -d rearr.xyz -d www.rearr.xyz -d gallery.rearr.xyz -d workflow.rearr.xyz -d flower.rearr.xyz -d shiny.rearr.xyz --nginx -m ljw2017@sjtu.edu.cn -n --cert-name wulab --test-cert

docker run -it --rm -v ./certbot/letsencrypt/:/etc/letsencrypt/ -v ./certbot/lib/letsencrypt/:/var/lib/letsencrypt/ certbot/certbot certonly -d qiangwulab.sjtu.edu.cn -d *.qiangwulab.sjtu.edu.cn -d rearr.xyz -d *.rearr.xyz --nginx -m ljw2017@sjtu.edu.cn -n --cert-name wulab --test-cert

# if ! grep -E "0 12 \* \* \* docker run -it --rm -v ./certbot/letsencrypt/:/etc/letsencrypt/ -v ./certbot/lib/letsencrypt/:/var/lib/letsencrypt/ certbot/certbot renew -q" /etc/crontab
# then
#     echo "0 12 * * * docker run -it --rm -v ./certbot/letsencrypt/:/etc/letsencrypt/ -v ./certbot/lib/letsencrypt/:/var/lib/letsencrypt/ certbot/certbot renew -q" >>/etc/crontab
# fi
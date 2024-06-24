#!/bin/bash

dhcpIp=$(ip a show dev enp5s0 | grep "inet " | awk '{print $2}' | cut -d'/' -f1)
printf "[Service]\nEnvironment=\"HTTP_PROXY=socks5://%s:1081\"\nEnvironment=\"HTTPS_PROXY=socks5://%s:1081\"" $dhcpIp $dhcpIp >$HOME/.config/systemd/user/docker.service.d/http-proxy.conf
systemctl --user daemon-reload
systemctl --user restart docker
#!/bin/bash

dirpath='/etc/systemd/system/docker.service.d/'
filename='http-proxy.conf'
sudo mkdir -p $dirpath
port='[Service]Environment="HTTP_PROXY=http://proxy.example.com:80/"'
echo $port > $dirpath$filename
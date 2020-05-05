FROM ubuntu:18.04

RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections

RUN apt-get -y update && apt-get -y install default-jdk wget bcftools

RUN wget https://github.com/broadinstitute/picard/releases/download/2.22.4/picard.jar \
  -P /opt/
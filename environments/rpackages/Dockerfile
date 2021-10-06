FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get clean && apt-get update --fix-missing && apt-get -y install unzip r-base curl libxml2-dev libcurl4-openssl-dev libssl-dev

RUN R -e 'install.packages(c("cowplot", "plotly", "htmlwidgets", "qqman", "tidyverse"))'

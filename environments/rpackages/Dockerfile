FROM rocker/tidyverse:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y install unzip curl libxml2-dev libcurl4-openssl-dev libssl-dev pandoc

RUN R -e 'install.packages(c("cowplot", "plotly", "htmlwidgets", "qqman"))'

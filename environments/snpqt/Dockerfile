FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get clean && apt-get update --fix-missing && apt-get -y install unzip r-base curl libxml2-dev libcurl4-openssl-dev libssl-dev

RUN curl https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 \
  -L -o /opt/bcftools.tar.bz2
RUN cd /opt/ && tar -xjf bcftools.tar.bz2
RUN cd /opt/bcftools-1.10.2/ && ./configure && make && make install

ENV BCFTOOLS_PLUGINS=/opt/bcftools-1.10.2/plugins

RUN R -e 'install.packages(c("cowplot", "plotly", "htmlwidgets", "qqman", "tidyverse"))'

RUN curl https://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip \
  --output /opt/plink2.zip
RUN unzip /opt/plink2.zip -d /opt/plink2/
ENV PATH="/opt/plink2:${PATH}"

RUN curl http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip \
  --output /opt/plink.zip
RUN unzip /opt/plink -d /opt/plink
ENV PATH="/opt/plink:${PATH}"

RUN apt-get update --fix-missing && apt-get -y install samtools tabix pandoc

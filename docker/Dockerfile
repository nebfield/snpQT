FROM ubuntu:18.04

RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections

RUN apt-get -y update && apt-get -y upgrade && apt-get -y install curl unzip \
  python-pip zlib1g-dev libbz2-dev liblzma-dev tabix r-base r-base-dev \
  libcurl4-openssl-dev libssl-dev libxml2-dev

RUN curl http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200219.zip \
  --output /opt/plink.zip
RUN unzip /opt/plink.zip -d /opt/plink
ENV PATH="/opt/plink:${PATH}"

RUN pip install snpflip biopython

RUN curl https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 \
  -L -o /opt/bcftools.tar.bz2
RUN cd /opt/ && tar -xf bcftools.tar.bz2
RUN cd /opt/bcftools-1.10.2/ && ./configure && make && make install

ENV BCFTOOLS_PLUGINS=/opt/bcftools-1.10.2/plugins

RUN R -e 'install.packages("tidyverse")'
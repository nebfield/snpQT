FROM ubuntu:18.04

RUN apt-get -y update && apt-get -y upgrade && apt-get -y install curl unzip 

RUN curl http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20191130.zip --output /opt/plink.zip 
RUN unzip /opt/plink.zip -d /opt/plink
ENV PATH="/opt/plink:${PATH}"

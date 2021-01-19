FROM ubuntu:18.04

RUN apt-get update && apt-get -y install unzip 
ADD impute_v5.1.zip /opt/
RUN ls /opt/ 
RUN unzip /opt/impute_v5.1.zip -d /opt/
ENV PATH="/opt/impute_v5.1:${PATH}"
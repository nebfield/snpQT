FROM ubuntu:18.04

RUN apt-get update && apt-get -y install unzip 
ADD impute5_v1.1.4.zip /opt/
RUN ls /opt/ 
RUN unzip /opt/impute5_v1.1.4.zip -d /opt/
RUN mkdir -p /opt/impute5
RUN mv /opt/impute5_v1.1.4/impute5_1.1.4_static /opt/impute5/impute5
RUN mv /opt/impute5_v1.1.4/imp5Converter_1.1.4_static /opt/impute5/imp5Converter
RUN chmod +x /opt/impute5/imp5Converter /opt/impute5/impute5
ENV PATH="/opt/impute5/:${PATH}"

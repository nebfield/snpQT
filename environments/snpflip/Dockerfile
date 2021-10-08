FROM ubuntu:xenial-20200619

RUN apt-get update && apt-get -y install python curl unzip
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py && python get-pip.py 
RUN pip2 install pandas==0.21.0 numpy==1.16.6 snpflip

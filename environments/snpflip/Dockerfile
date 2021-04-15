FROM ubuntu:xenial-20200619

RUN apt-get update && apt-get -y install python curl unzip
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py && python get-pip.py 
RUN pip2 install pandas==0.21.0 snpflip
RUN curl http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200616.zip --output /opt/plink.zip 
RUN unzip /opt/plink.zip -d /opt/plink/
ENV PATH="/opt/plink:${PATH}"


# https://github.com/odelaneau/shapeit4/issues/3#issuecomment-509376555 
FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive
ARG SAMTOOLS_VERSION=1.9
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ca-certificates \
                 libbz2-dev \
                 libcurl4-openssl-dev \
                 liblzma-dev \
                 libncurses5-dev \
                 autoconf \
                 automake \
                 bzip2 \
                 gcc \
		 g++ \
                 make \
                 wget \
                 zlib1g-dev \
		 libboost-all-dev \
		 git && \
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && ./configure CFLAGS="-static" --without-curses && make -s all all-htslib && make install install-htslib && cd - && \ 
    rm -rf samtools-${SAMTOOLS_VERSION}* && \
    git clone https://github.com/odelaneau/shapeit4.git && \
    cd shapeit4 && \
    cat makefile | sed 's|HTSLIB_INC=$(HOME)/Tools/htslib-1.9|HTSLIB_INC=/usr/local/include|' | sed 's|HTSLIB_LIB=$(HOME)/Tools/htslib-1.9/libhts.a|HTSLIB_LIB=/usr/local/lib/libhts.a|' > makefile2 && \
    mv makefile makefile.old && \
    mv makefile2 makefile && \
    make

ENV PATH="/shapeit4/bin/:${PATH}"

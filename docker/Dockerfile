FROM ubuntu:18.04
MAINTAINER Trevor Pesout tpesout@ucsc.edu

ARG git_commit

# update and install dependencies
RUN apt-get update && \
    apt-get -y install time git make wget autoconf gcc g++ zlib1g-dev libcurl4-openssl-dev libbz2-dev libhdf5-dev liblzma-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# install cmake
WORKDIR /tmp
RUN mkdir /opt/cmake && \
    wget https://cmake.org/files/v3.11/cmake-3.11.4-Linux-x86_64.sh && \
    sh /tmp/cmake-3.11.4-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && \
    ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake

# get samtools
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar xvf samtools-1.9.tar.bz2 && \
    rm -r /opt/samtools/samtools-1.9.tar.bz2 && \
    cd samtools-1.9/ && \
    autoheader && \
    autoconf -Wno-header && \
    ./configure --without-curses --disable-lzma && \
    make && \
    ln -s /opt/samtools/samtools-1.9/samtools /usr/local/bin/samtools

# get marginPolish
WORKDIR /opt
RUN git clone https://github.com/UCSC-nanopore-cgl/margin.git && \
    cd /opt/margin && \
    git fetch && \
    git checkout $git_commit && \
    git submodule update --init && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    rm -rf /opt/margin/.git && \
    ln -s /opt/margin/build/margin /usr/local/bin/margin && \
    margin phase -h && \
    margin polish -h

WORKDIR /data

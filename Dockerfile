# Get the base Ubuntu image from Docker Hub
FROM ubuntu:20.04

# Install GCC and dependencies
RUN apt-get -y update \
    && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata \
    && apt-get install -y git wget bzip2 autoconf automake make cmake gcc g++ perl zlib1g-dev libbz2-dev liblzma-dev \
    libcurl4-gnutls-dev libssl-dev libncurses5-dev

WORKDIR /usr/src/
RUN wget https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2 -O htslib.tar.bz2 \
    && tar -xjvf htslib.tar.bz2 \
    && cd htslib-1.15.1 \
    && autoreconf -i \
    && ./configure \
    && make \
    && make install

# Start building
COPY . /usr/src/mvtool
WORKDIR /usr/src/mvtool/build
RUN cmake -DENABLE_MIMALLOC=ON .. \
    && make -j

# Get binaries
WORKDIR /mvtool/bin
RUN cp \
    /usr/src/mvtool/build/mvtool \
    .
ENV PATH /mvtool/bin:$PATH
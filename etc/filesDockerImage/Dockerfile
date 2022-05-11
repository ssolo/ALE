# Dockerfile for ALE with BPP v2.4.1
# Version: 
#	- v0.1 (2022-05-01)

FROM ubuntu:22.04

RUN apt-get update && \
    apt-get clean && \
    apt-get install -qy \
			cmake \
			libboost-all-dev \
			g++-11 \
			git \
			make \
			python3 \
			wget \
			ca-certificates \
			openssl \
			build-essential \
			libeigen3-dev


WORKDIR /opt

# Install libboost
RUN mkdir -p {bpp/bpp-core-build,bpp/bpp-seq-build,bpp/bpp-phyl-build}

# Install Bio++ v2.4.1
WORKDIR /opt/bpp

# pull repositories
RUN git clone https://github.com/BioPP/bpp-core
RUN git clone https://github.com/BioPP/bpp-seq
RUN git clone https://github.com/BioPP/bpp-phyl

# freeze to version 2.4.1
WORKDIR /opt/bpp/bpp-core
RUN git checkout tags/v2.4.1 -b version2.4.1
RUN sed -i '45i#include <limits>' src/Bpp/Graph/GlobalGraph.cpp

WORKDIR /opt/bpp/bpp-seq
RUN git checkout tags/v2.4.1 -b version2.4.1

WORKDIR /opt/bpp/bpp-phyl
RUN git checkout tags/v2.4.1 -b version2.4.1

WORKDIR /opt/bpp/bpp-core-build
RUN cmake ../bpp-core -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_TESTING=FALSE
RUN make
RUN make install

WORKDIR /opt/bpp/bpp-seq-build
RUN cmake ../bpp-seq -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_TESTING=FALSE
RUN make
RUN make install

WORKDIR /opt/bpp/bpp-phyl-build
RUN cmake ../bpp-phyl -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_TESTING=FALSE
RUN make
RUN make install


# Install ALE
WORKDIR /opt
RUN mkdir ALE-build

RUN ln -s /usr/include/eigen3/Eigen /usr/include/Eigen

RUN git clone https://github.com/ssolo/ALE.git

WORKDIR /opt/ALE
RUN git checkout bppv241

WORKDIR /opt/ALE-build
RUN cmake ../ALE -DCMAKE_LIBRARY_PATH=/usr/local -DCMAKE_INCLUDE_PATH=/usr/local
RUN make

RUN for binary in $PWD/bin/*; do ln -s $binary /usr/local/bin/; done

ENV LD_LIBRARY_PATH /usr/local/lib

FROM debian:jessie

MAINTAINER Bastien Boussau bastien.boussau@univ-lyon1.fr

RUN apt-get update && \
    apt-get clean && \
    apt-get install --no-install-recommends -qy \
                      cmake \
                      g++-4.9 \
                      git \
                      libboost-all-dev \
                      make \
                      python3 \
                      wget

#### install bpp
WORKDIR /opt
RUN echo "deb http://biopp.univ-montp2.fr/repos/apt/ Trusty main" >> /etc/apt/sources.list; \
    wget http://biopp.univ-montp2.fr/repos/apt/conf/biopp.gpg.key && \
    apt-key add biopp.gpg.key && \
	  apt-get update && \
    apt-get install -qy libbpp-phyl-dev


RUN git clone git://github.com/ssolo/ALE /usr/local/ALE && \
    mkdir /usr/local/ALE/build

WORKDIR /usr/local/ALE/build

ENV LD_LIBRARY_PATH /usr/local/lib/

RUN cmake ..  -DCMAKE_CXX_COMPILER=/usr/bin/g++-4.9 && \
    make -j 4 && \
    for binary in $PWD/bin/*; do ln -s $binary /usr/local/bin/; done

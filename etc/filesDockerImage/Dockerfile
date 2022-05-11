# Dockerfile for ALE with BPP v3
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

# Install Bio++ v3 developmental version
WORKDIR /opt/bpp

# pull repositories
RUN git clone https://github.com/BioPP/bpp-core
RUN git clone https://github.com/BioPP/bpp-seq
RUN git clone https://github.com/BioPP/bpp-phyl

# freeze to version 3.0
WORKDIR /opt/bpp/bpp-core
RUN git checkout ded2725612e1bd850b1257ed11bd020fc21fc429

WORKDIR /opt/bpp/bpp-seq
RUN git checkout 577e1eb6392a9c0d69f571f423d182c5d82185c3

WORKDIR /opt/bpp/bpp-phyl
RUN git checkout eac896ed37c85fb4737b04f898c7a99c55662d1a
RUN sed -i "s/constexpr ExtendedFloat (const ExtendedFloat\& ef) noexcept : f\_ (ef.float_part()), exp\_ (ef.exponent_part())/ExtendedFloat (const ExtendedFloat\& ef) noexcept : f\_ (ef.f\_), exp\_ (ef.exp\_)/"  src/Bpp/Phyl/NewLikelihood/DataFlow/ExtendedFloat.h

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
RUN git checkout bppv3
RUN git pull

WORKDIR /opt/ALE-build
RUN cmake ../ALE -DCMAKE_LIBRARY_PATH=/usr/local -DCMAKE_INCLUDE_PATH=/usr/local -DBUILD_STATIC=OFF
RUN make

RUN for binary in $PWD/bin/*; do ln -s $binary /usr/local/bin/; done

WORKDIR /opt

ENV LD_LIBRARY_PATH /usr/local/lib


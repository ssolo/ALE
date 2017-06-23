#FROM ubuntu:16.04
FROM debian:jessie

MAINTAINER Bastien Boussau bastien.boussau@univ-lyon1.fr

RUN apt-get clean && apt-get update

RUN apt-get install --no-install-recommends -qy make \
    	    	    	wget \
			git \
			python3 \
			#-setuptools \
			#python-pip \


# to install phyldog
                        cmake \
                        g++-4.9 \
                        libboost-all-dev



#### install bpp
WORKDIR /opt
RUN echo "deb http://biopp.univ-montp2.fr/repos/apt/ Trusty main" >> /etc/apt/sources.list;\
    wget http://biopp.univ-montp2.fr/repos/apt/conf/biopp.gpg.key &&\
    apt-key add biopp.gpg.key && \
	apt-get update && \
    apt-get install -qy libbpp-phyl-dev


RUN git clone git://github.com/ssolo/ALE /usr/local/ALE &&  mkdir /usr/local/ALE/build
    
WORKDIR /usr/local/ALE/build

ENV LD_LIBRARY_PATH /usr/local/lib/
RUN cmake ..  -DCMAKE_CXX_COMPILER=/usr/bin/g++-4.9 && make -j 4
#########

#### build a working directory
#RUN mkdir -p /app
WORKDIR /app

#RUN ls /usr/local/ALE/build/bin

#RUN ls /usr/local/ALE/build

#RUN /usr/local/ALE/build/bin/ALEobserve

# Copy the current directory contents into the container at /app
ADD . /app


ENTRYPOINT ["python3.4", "launcher_script.py"]

#ENTRYPOINT ["/usr/local/ALE/build/bin/ALEobserve"]
#ENTRYPOINT ["/bin/cat"]
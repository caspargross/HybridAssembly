################################
#                              #
#         DOCKERFILE           #
#   Hybrid Assembly Pipeline   #
#                              #
################################

FROM ubuntu:16.04
MAINTAINER mail@caspar.one

# Commont applications
RUN apt-get update -y
RUN apt-get install -y python
RUN apt-get install -y perl

# Install JAVA
RUN add-apt-repository ppa:webupd8team/java &&\
    apt-get update &&\
    (echo "oracle-java8-installer shared/accepted-oracle-license-v1-1 select true" | sudo debconf-set-selections) &&\
    apt-get isntall -y oracle-java8-installer

# Install SPAdes


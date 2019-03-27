################################
#                              #
#         DOCKERFILE           #
#   Hybrid Assembly Pipeline   #
#                              #
################################
FROM continuumio/miniconda3
MAINTAINER Caspar Gross <mail@caspar.one>
LABEL description="contains all the dependencies for hybridAssembly pipeline at github.com/caspargross/hybridAssembly" 

# Install Java.
RUN \
  apt-get update && \
  apt-get install -y openjdk-8-jre gawk bc procps && \
  rm -rf /var/lib/apt/lists/*

# Set standard shell to bash
SHELL ["/bin/bash", "-ce"]

# Define commonly used JAVA_HOME variable
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64

# Install conda envs
ADD envs/ha_py36.yml /tmp/ha_py36.yml
RUN conda env create -f /tmp/ha_py36.yml -q && conda clean -a

# Fix Bug https://github.com/nf-core/ampliseq/issues/25 which prevents Bandage from drawing graphs because of missing qt lib
RUN echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

ADD envs/ha_py27.yml /tmp/ha_py27.yml
RUN conda env create -f /tmp/ha_py27.yml -q && conda clean -a


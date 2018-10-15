################################
#                              #
#         DOCKERFILE           #
#   Hybrid Assembly Pipeline   #
#                              #
################################
FROM continuumio/miniconda
MAINTAINER Caspar Gross <mail@caspar.one>
LABEL description="contains all the dependencies for hybridAssembly pipeline at github.com/caspargross/hybridAssembly" 

# Install basic packages into docker container
RUN apt-get update && apt-get install procps bc gawk 

# Install java
# Install Java.
RUN \
  echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | debconf-set-selections && \
  add-apt-repository -y ppa:webupd8team/java && \
  apt-get update && \
  apt-get install -y oracle-java8-installer && \
  rm -rf /var/lib/apt/lists/* && \
  rm -rf /var/cache/oracle-jdk8-installer

# Define working directory.
WORKDIR /data

# Define commonly used JAVA_HOME variable
ENV JAVA_HOME /usr/lib/jvm/java-8-oracle

# Set standard shell to bash
SHELL ["/bin/bash", "-c"]

# Install conda environments
COPY envs/ha_py36.yml /
RUN conda env create -f /ha_py36.yml -q && conda clean -a

COPY envs/ha_py27.yml /
RUN conda env create -f /ha_py27.yml -q && conda clean -a

# Download checkM database
RUN mkdir -p /data && cd /data
RUN wget -q -O checkm_data.tar.gz https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz &&  tar -zxvf checkm_data.tar.gz && rm checkm_data.tar.gz && source activate ha_py27 && checkm_db="/data" &&  echo ${checkm_db} | checkm data setRoot ${checkm_db}

# Download CARD-Antibiotic resistance database
RUN wget -q -O card-data.tar.bz2 https://card.mcmaster.ca/latest/data && tar xfvj card-data.tar.bz2 && bash && source activate ha_py27 && rgi load --afile card.json


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

# Define working directory.
# WORKDIR /data

# Install conda envs
ADD envs/ha_py36.yml /tmp/ha_py36.yml
RUN conda env create -f /tmp/ha_py36.yml -q && conda clean -a

ADD envs/ha_py27.yml /tmp/ha_py27.yml
RUN conda env create -f /tmp/ha_py27.yml -q && conda clean -a

# Download checkM database
RUN mkdir -p /data && cd /data
RUN wget -q -O checkm_data.tar.gz https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz &&  tar -zxvf checkm_data.tar.gz && rm checkm_data.tar.gz && source activate ha_py27 && checkm_db="/data" &&  echo ${checkm_db} | checkm data setRoot ${checkm_db}

# Download CARD-Antibiotic resistance database
RUN wget -q -O card-data.tar.bz2 https://card.mcmaster.ca/latest/data && tar xfvj card-data.tar.bz2 && bash && source activate ha_py27 && rgi load --afile card.json

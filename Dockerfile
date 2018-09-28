################################
#                              #
#         DOCKERFILE           #
#   Hybrid Assembly Pipeline   #
#                              #
################################

FROM continuumio/miniconda
MAINTAINER Caspar Gross <mail@caspar.one>
LABEL description="contains all the dependencies for hybridAssembly pipeline at github.com/caspargross/hybridAssembly" 


# Install conda environments
COPY envs/ha_py36.yml /
RUN conda env create -f /ha_py36.yml -q

COPY envs/ha_py27.yml /
RUN conda env create -f /ha_py27.yml -q

RUN conda clean -a
RUN bash

# Download checkM database
RUN mkdir -p /data && cd /data
RUN source activate ha_py27 
RUN wget -O checkm_data.tar.gz https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
RUN tar -zxvf checkm_data.tar.gz 
RUN rm checkm_data.tar.gz
RUN printf '/data' | checkm data setRoot

# Download CARD-Antibiotic resistance database
RUN wget -O card-data.tar.bz2 https://card.mcmaster.ca/latest/data
RUN tar xfvj card-data.tar.bz2
RUN rgi load --afile card.json


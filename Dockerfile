################################
#                              #
#         DOCKERFILE           #
#   Hybrid Assembly Pipeline   #
#                              #
################################

FROM continuumio/miniconda
MAINTAINER Caspar Gross <mail@caspar.one>
LABEL description="contains all the dependencies for hybridAssembly pipeline at github.com/caspargross/hybridAssembly" 

# Set standard shell to bash
SHELL ["/bin/bash", "-c"]

# Install conda environments
COPY envs/ha_py36.yml /
RUN conda env create -f /ha_py36.yml -q && conda clean -a

COPY envs/ha_py27.yml /
RUN conda env create -f /ha_py27.yml -q && conda clean -a

# Download checkM database
RUN mkdir -p /data && cd /data
RUN wget -q -O checkm_data.tar.gz https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz &&  tar -zxvf checkm_data.tar.gz && rm checkm_data.tar.gz && source activate ha_py27 && printf '/data\n' | checkm data setRoot

# Download CARD-Antibiotic resistance database
RUN wget -q -O card-data.tar.bz2 https://card.mcmaster.ca/latest/data && tar xfvj card-data.tar.bz2 && bash && source activate ha_py27 && rgi load --afile card.json


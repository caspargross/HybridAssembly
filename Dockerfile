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


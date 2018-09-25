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
RUN conda env create -f /ha_py36.yml && conda clean -a

COPY envs/ha_py27.yml /
RUN conda env create -f /ha_py27.yml && conda clean -a

COPY envs/ha_pl.yml /
RUN conda env create -f /ha_pl.yml && conda clean -a

COPY envs/ha_quast.yml /
RUN conda env create -f /ha_quast.yml && conda clean -a

COPY envs/ha_uni.yml /
RUN conda env create -f /ha_uni.yml && conda clean -a


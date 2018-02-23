hybridAssembly
===============

hybridAssembly is a pipeline for hybrid de novo assembly of bacterial genomes. The pipeline is based around a [NextFlow](https://github.com/nextflow-io/nextflow) script. 

Aims of this pipeline:
* Give an overview about different approaches to hybrid Assembly
* Compare the results of different tools
* Learning worklfow management with Nextflow

This pipeline is mainly intended for personal use and therefore not well documented or easy to implement. If I find the time, I might create a docker container to share all tools used. 

Installation:
-------------

1) Clone this repository

2) Install Nextflow [NextFlow](https://github.com/nextflow-io/nextflow)

3) Change the path variables for the different tools in the config file  'nextflow.config'

Assembly paths:
---------------

<img src="./FlowSchema.svg">

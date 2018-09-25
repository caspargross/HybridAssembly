hybridAssembly
===============

hybridAssembly is a pipeline for hybrid de novo assembly of bacterial genomes. The pipeline uses a [NextFlow](https://github.com/nextflow-io/nextflow) script to combine several tools. 

Aims of this pipeline:
* Give an overview about different approaches to hybrid Assembly
* Compare the results of different tools
* Using Nextflow as a workflow management tool.

Currently this pipeline is mainly intended for personal use and therefore not well documented and in final state. Hopefully there will be a docker image including all the dependencies soon.  

Installation:
-------------

1) Clone this repository

2) Install Nextflow [NextFlow](https://github.com/nextflow-io/nextflow)

3) Change the path variables for the different tools in the config file  'nextflow.config'


Dependencies:
-------------

For all needed tools see the environment variables in the config file. For some tools it is recommended to install conda packages.


Assembly paths:
---------------
There are seven different run configuration which select different assembly tools:

* `spades_sspace` SPAdes assembly with SSPACE scaffolding
* `spades_links` SPAdes assembly with LINKS scaffolding
* `canu` Canu assembly
* `unicycler` [Unicycler](https://github.com/rrwick/Unicycler) full bacterial assembly pipeline
* `flye` [Flye](https://github.com/fenderglass/flye) assembler
* `miniasm` [Miniasm](https://github.com/lh3/miniasm) assembler
* `all` Execute all paths in parallel to compare results

<img src="./FlowSchema.svg">

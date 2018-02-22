HybridAssembly
===============

Nextflow pipeline for hybrid assembly of bacterial genomes using a combination of long and short read sequencing data.

Workflow:
    1) *Assembly* using SPAdes hybrid assembly
    2) *Scaffolding*: Although SPAdes has his own scaffolding implementation, the results are optimized with either SSPACE or LINKS.
    3) Gapfilling: 
    4) Alignment and variant calling using Mummer.


![FlowDiagram](https://raw.githubusercontent.com/caspargross/hybridAssembly/master/FlowSchema.svg)

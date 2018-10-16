Requirements
------------

### Operating System

Currently nextflow is only availble for POSIX compliant systems (Linux, MacOS ...)

### Software

- Java version > 8.0
- Nextflow, current version
- Docker (recommended)
- Singularity, if Docker is not possible
- Conda (optional) if Docker/Singularity is not an option

see Installation instructions below

### Hardware

This pipeline has only been tested on a Server it is recommended to have at leas 32Gb of memory and > 8 CPU. The most demanding step (checkM) can be skipped by chosing `--fast` as a run parameter

The test profile (start with options `-profile test`) runs fine on a personal Laptop

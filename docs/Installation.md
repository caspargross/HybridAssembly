[Nextflow](https://www.nextflow.io/) needs to be installed to run this pipeline. Nextflow can be used on any POSIX compatible system (Linux, OS X, etc). It requires BASH and Java 8 (or higher) to be installed.

Nexflow Installation:
---------------------
Execute the following steps to install nextflow on your machine

```
# Check java version (needs 1.8 or higher)
java -version

# Download and install nextflow
wget -qO- get.nextflow.io | bash

# (Optional) move the file to a $PATH accessible directory
sudo mv nextflow /usr/local/bin/
```

There are several ways to install the Tools needed by the pipeline:

Docker installation:
--------------------

Using the docker container system is the most convenient way to run this pipeline.  Install docker using the package manager of your system or look on their (homepage)[https://docs.docker.com/install/] for more specific installation procedures.

Singularity installation:
-------------------------

If runnign docker on your machine is not possible or desirable (i.e. on high performance cluster machines) you can use [Singularity](https://singularity.lbl.gov/install-request nstead). See the homepage for installation procedures.

Local installation using conda:
-------------------------------

These steps are only necessary when you dont want to use docker/singularity. All the tools needed are contained in two [conda](https://conda.io/docs/) environments located in the `env/` Folder. Install them so you can call them from you path later.

1) Install miniconda (skip if already available). Instructions [here](https://conda.io/docs/user-guide/install/index.html#)

2) Install environments
    
``` 
conda create -f envs/ha_py27.yml
conda create -f envs/ha_py36.yml

# Optional (removes temp files)
conda clean -a
```

3) Download required databases

Download and install [checkM](https://github.com/Ecogenomics/CheckM) Database:

```
mkdir -p /checkm_data && cd /checkm_data
wget -q -O checkm_data.tar.gz https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz 
tar -zxvf checkm_data.tar.gz && rm checkm_data.tar.gz

source activate ha_py27 
echo "checkm_data" | checkm data setRoot "checkm_data"
source deactivate ha_py27
```

Download and include CARD-Antibiotic resistance database

```
wget -q -O card-data.tar.bz2 https://card.mcmaster.ca/latest/data
tar xfvj card-data.tar.bz2

source activate ha_py27 
rgi load --afile card.json 
source deactivate ha_py27
```


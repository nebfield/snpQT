# Installation

[TOC]

`snpQT` uses docker containers, Anaconda and/or environment modules to install and run an underlying collection of bioinformatics software. Before you can download and run the pipeline you'll need:

* [Nextflow](https://www.nextflow.io) `>=20.10.0`
* [Docker](https://docs.docker.com/get-docker/), [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) or [Environment Modules](http://modules.sourceforge.net/)
    * If using docker, you'll need to [run the post-install steps](https://docs.docker.com/engine/install/linux-postinstall/)

!!!Tip
	* Conda is suitable for users who are not interested in performing local imputation and who do not have root access in their machines. However, with Conda you can still run imputation-related workflows, like pre-imputation and post-imputation QC, as well as all the rest QC-related workflows of snpQT.
	* Docker requires root access, while enables the installation of `impute5`, the imputation software.
	* Environment modules are useful to run in HPC environments, where root access is not available, however they still require some user configuration because installed packages and package names are custom to each HPC environment
	
* A reasonably powerful computer running something Linux flavoured 
    * `snpQT` has been tested on centOS and Ubuntu

If you can run `nextflow run hello`, `docker run hello-world` (or `conda -V`) without any error messages you should be ready to download `snpQT`. 
    
!!! Warning
	* In order to make sure you have the latest version of nextflow you can run `nextflow self-update`
	* Docker requires root (superuser / administrative) access to run. If you don't have root access you will need to run `snpQT` using conda. With conda, `snpQT` can't do imputation.

	
## Install via source

To begin the installation, clone the repository:

```
git clone https://github.com/nebfield/snpQT.git
cd snpQT
nextflow run main.nf
```

You should see some helpful usage information printed to your terminal if everything went well.

## Core reference data

To do anything interesting with `snpQT` you'll need to download some reference data. The core reference data is required for build conversion, sample and variant quality control, population stratification, pre-imputation and post-imputation.

```
nextflow run main.nf -profile [conda/docker] --download_db core
```

This will download the core reference files, prepare them and put them in a `./db/` folder in the `snpQT` directory. On our computers it takes around an hour to run but this may take longer depending on your network. The core dataset requires about ~43GB of initial storage of intermediate files in the `./work/` that can be removed using `nextflow run main.nf --download_db core && rm -r work` and only 19.7Gb of reference files that are stored in database directory `./db/`.

Another way to download the required reference data is to directly download an already processed `.tar.gz` file (17.3Gb) from one of our servers, which on our computers takes around 40 minutes to finish, using the following lines of code:

```
mkdir ./db/
cd db
wget https://sys-myo.com/snpQT/core.tar.gz
tar -xvf core.tar.gz
```

!!!Tip
	If you are not interested in using pre-imputation (and local imputation see below), you can significantly free up space by deleting `All_20180423.vcf.gz` (and `All_20180423.vcf.gz.tbi`), the largest file in the database (14.3Gb) which is only used for this purpose and it is not needed in the other workflows. This will leave you with only ~5Gb size storage of reference files.

## Imputation reference data

To do local imputation with `snpQT` you'll need to download additional reference data after downloading the core reference data. 

```
nextflow run main.nf -profile [conda/docker] --download_db impute
```

This will download the imputation reference files and put them in `snpQT/db/impute`. Downloading and processing the imputation data can take an hour or two depending on your network and computer. These reference files require an additional 30GB of storage.

Another way to download the required reference data is to directly download an already processed `.tar.gz` file (13Gb) from one of our servers, in our computers the download takes about 25 minutes to finish. For this purpose, you can use the following lines of code:

```
mkdir ./db/
cd db
wget https://sys-myo.com/snpQT/impute.tar.gz
tar -xvf impute.tar.gz
```

## The `work` directory

Pipelines like `snpQT` have many steps, and make a lot of files that aren't shown directly to the user. Nextflow stores all of these intermediate files in the `snpQT/work/` directory by default ([this can be changed with `nextflow -w`](https://www.nextflow.io/docs/latest/cli.html)). The work directory can get **very big very quickly**. After downloading your data you might want to delete this folder, and regularly check the size of this folder over time.

`snpQT` will remember previous work it's done. This means if you want to run a different stage of the pipeline on the same input data, it should skip work already done. If you delete the `work/` folder you will lose this cache. 

!!! warning
    If you delete the `work` directory you will lose:

      * Conda environments that have been automatically set up
      * Any cached work
      
    You might want to try running [`nextflow clean`](https://www.nextflow.io/docs/latest/cli.html) instead 

## Imputation prerequisite 

`snpQT` uses Docker biocontainers or conda to provide bioinformatics software for each part of the pipeline. However, `impute5` has a restrictive license and is only free for academic users:

<https://jmarchini.org/software/#impute-5>

If you want to run imputation locally on your computer, you will need:

* A copy of `impute5` (version 5.1) in the `environments/impute5/` directory
* To manually build the imputation docker image

```
cd environments/impute5/
docker build -t snpqt/imputation .
```

When the Dockerfile builds successfully you'll be able to do local imputation.

## Run the toy dataset

Once installed, it's sensible to test `snpQT` with the provided toy dataset. Click Next to go to the Quickstart section. A more detailed [Tutorial](https://tutorial-snpqt.readthedocs.io/en/latest/user-guide/results/) section is provided later.

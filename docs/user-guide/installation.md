# Installation

[TOC]

## Requirements and dependencies

The mandatory software requirements to run `snpQT` are quite simple:

* [Nextflow](https://www.nextflow.io) `>=20.10.0` 
* A computer running something modern and Linux flavoured (we've tested
  on CentOS)

If you can run `nextflow run hello` without any error messages you should be
ready to download `snpQT`.

The hardware requirements depend on the size of your input data and analysis
task. Quality control on hundreds or thousands of samples can be comfortably
done on a laptop. Imputation of thousands of samples will require a more
powerful server or high-performance computing cluster. We've measured RAM usage
of up to 80GB RAM per chromosome on a dataset of 40,000 individuals when doing
imputation. Reference data can take up to 50GB of disk space.

`snpQT` would run OK if all of the underlying bioinformatics software we use
happened to be already correctly installed and configured on your computer. This
is unlikely, so we've provided a few methods, which nextflow calls `-profiles`,
to try and make this as easy as possible. It's very important to install:

* [Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

And at least one of these, depending on your local setup, with Singularity being
the best option:

* [Singularity](https://sylabs.io/guides/3.3/user-guide/index.html) 
* [Docker](https://docs.docker.com/get-docker/)
* [Environment Modules](http://modules.sourceforge.net/)

!!! Tip 
    To do imputation locally you'll need to use either singularity, docker, or
    environment modules. All of these methods require some [additional
    configuration](#configuring-imputation-optional-but-a-good-idea).
	
## Install latest release

The [release page always contains the latest version of
`snpQT`](https://github.com/nebfield/snpQT/releases). For example, you might do
something like:

```
$ wget https://github.com/nebfield/snpQT/archive/refs/tags/v0.1.3.zip
$ unzip v0.1.3.zip
$ cd snpQT-0.1.3/
$ nextflow run main.nf
```

The exact syntax will change depending on what version of `snpQT` you've
downloaded. If everything went well you should see:

```
=================================================================
snpQT is ready to make your single-nucleotide polymorphisms cute!
v0.1.3 - Fluffy penguin, 2021-06-15
... (helpful stuff)
```

## Download pre-built reference data

To do anything interesting with `snpQT` you'll need to download some reference
data. The core reference data is required for build conversion, sample and
variant quality control, population stratification, pre-imputation and
post-imputation.

`snpQT` has two reference data sets depending on your analysis task at
[10.5281/zenodo.4916468](https://doi.org/10.5281/zenodo.4916468):

* `core.tar.gz` is about 18.6GB, and required for all `snpQT` workflows except
  imputation 
* `impute.tar.gz` is about 15.0GB and is required to do imputation

To set up the core data:

```
$ mkdir snpQT-0.1.3/db
$ wget 'https://zenodo.org/record/4916469/files/core.tar.gz?download=1' -O snpQT-0.1.3/db/core.tar.gz
$ cd snpQT-0.1.3/db && tar -xvf core.tar.gz
```

And optionally if you're doing imputation:

```
$ wget 'https://zenodo.org/record/4916469/files/impute.tar.gz?download=1' -O snpQT-0.1.3/db/impute.tar.gz
$ cd snpQT-0.1.3/db && tar -xvf impute.tar.gz 
```

## Build your own reference data (optional)

If you're feeling masochistic you can build your own reference data. This will
take a decent chunk of computing power and time:

```
nextflow run main.nf -profile conda --download_db core
```

This will download the core reference files, prepare them and put them in a
`db/` folder in the `snpQT` directory. On our computers it takes around an hour
to run but this may take longer depending on your network. The core dataset
requires about ~43GB of initial storage of intermediate files in the `work/`
that can be removed using `nextflow run main.nf --download_db core && rm -r
work` and only 19.7GB of reference files that are stored in database directory
`db/`.


### The `work` directory

Pipelines like `snpQT` have many steps, and make a lot of files that aren't
shown directly to the user. Nextflow stores all of these intermediate files in
the `snpQT/work/` directory by default (this can be changed with `nextflow -w`). The
work directory can get **very big very quickly**. After downloading your data
you might want to delete this folder, and regularly check the size of this
folder over time.

`snpQT` will remember previous work it's done. This means if you want to run a
different stage of the pipeline on the same input data, it should skip work
already done. If you delete the `work/` folder you will lose this cache.

!!! warning

    If you delete the `work` directory you will lose:

      * Conda environments that have been automatically set up
      * Any cached work
      
    You might want to try running [`nextflow clean`](https://www.nextflow.io/docs/latest/cli.html) instead

## Configuring imputation (optional, but a good idea)

`snpQT` uses Anaconda and containers to reproducibly install and configure
bioinformatics software for each part of the pipeline. However, `impute5` has a
restrictive license and is only free for academic users:

* [https://jmarchini.org/software/#impute-5](https://jmarchini.org/software/#impute-5)
* [https://jmarchini.org/licence/](https://jmarchini.org/licence/)

Imputation takes a lot of computing power, so is ideally suited for powerful
servers or HPCs. If you want to configure impute5 to run in a container,
you will need:

* A copy of `impute5` (version 5.1.1.4) in the `environments/impute5/` directory
* To manually build the imputation singularity definition file (see
  `environments/impute5/impute5.def`) or dockerfile (see
  `environments/impute5/Dockerfile`) on  a machine with root access

Assuming you're in the `snpQT` base directory:

```
$ cd environments/impute5/
$ singularity build impute5.sif impute5.def
```

!!! Tip

    You will need root permissions to build singularity containers. You can
    build the container on a machine you have root access on and transfer it
    (the `.sif` file) to your HPC, or use the very useful [remote
    builder](https://cloud.sylabs.io/builder) service.

You might want to configure impute5 using environment modules, which are often
installed on HPCs. However, software naming schemes can be
inconsistent across environment modules installations on different HPC
systems. You will probably need to modify the module parameter for imputation processes in
`conf/modules.config`. On our HPC, `module avail` lists impute5 as
`impute/5_1.1.4`, which is the default:

```
withName: 'convert_imp5|impute5' {
    module = 'impute/5_1.1.4'
  }
```

Changing this to whatever your installation of impute5 v5.1.1.4 is called will
let you run imputation locally. We recommend Singularity because containers are
cool!

## Run the toy dataset

Once installed, it's sensible to test `snpQT` with the provided toy
dataset. Click Next to go to the Quickstart section. A more detailed
[Tutorial](https://tutorial-snpqt.readthedocs.io/en/latest/user-guide/results/)
section is provided later.

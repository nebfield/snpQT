# Advanced installation

[TOC]

## Clone development version

Download the most recent development version (pulling all the commits in our repository):

```
$ git clone https://github.com/nebfield/snpQT.git
```

`main.nf` is the entrypoint of the pipeline (e.g. `nextflow run main.nf`).

The development version may break randomly.

## Additional profiles

In total there are four profiles that control how to run modules. We recommend:

* `-profile conda`
* `-profile singularity`

Which were described in the [profiles](../quickstart/profiles.md) section. 

We also provide:

* `-profile docker`
* `-profile modules`

[Docker](https://docs.docker.com/get-docker/) requires root (superuser) access
to build and run containers. We also assume that you have run the [post-install
steps for Linux](https://docs.docker.com/engine/install/linux-postinstall/). 

Modules requires Anaconda and [Environment
modules](https://modules.readthedocs.io/en/latest/) to be installed, and is
mostly useful for cluster environments if you can't run Singularity.

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

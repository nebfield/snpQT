# Advanced installation

[TOC]

## Clone development version

```
$ git clone https://github.com/nebfield/snpQT.git
```

`main.nf` is the entrypoint of the pipeline (e.g. `nextflow run main.nf`).

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

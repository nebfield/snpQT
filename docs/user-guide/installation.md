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


## The `work` directory

Pipelines like `snpQT` have many steps, and make a lot of files that aren't
shown directly to the user. Nextflow stores all of these intermediate files in
the `snpQT/work/` directory by default (this can be changed with `nextflow
-w`). The work directory can get **very big very quickly**. After downloading
your data you might want to delete this folder, and regularly check the size of
this folder over time.

`snpQT` will remember previous work it's done. This means if you want to run a
different stage of the pipeline on the same input data, it should skip work
already done. If you delete the `work/` folder you will lose this cache.

!!! warning

    If you delete the `work` directory you will lose:

      * Conda environments that have been automatically set up
      * Any cached work
      
    You might want to try running [`nextflow clean`](https://www.nextflow.io/docs/latest/cli.html) instead



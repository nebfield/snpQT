# Installation

Assuming you're running a modern version of Linux, firstly install
[Nextflow](https://www.nextflow.io):

```
$ curl -s https://get.nextflow.io | bash 
```

Then simply run:

```
$ nextflow run nebfield/snpqt -r v0.1.4
```

Where `-r v0.1.4` reflects the latest release of `snpQT`. You can find different
releases of snpQT on our [releases
page](https://github.com/nebfield/snpQT/releases).  You should see a helpful
message in your terminal if everything went well.

Finally, you must download some reference data
([10.5281/zenodo.4916468](https://doi.org/10.5281/zenodo.4916468)) for snpQT to
work:

```
$ nextflow info nebfield/snpqt
```

Download the reference files to the local path (by default
`$HOME/.nextflow/assets/nebfield/snpqt`):

```
$ cd $HOME/.nextflow/assets/nebfield/snpqt
$ mkdir db
$ wget 'https://zenodo.org/record/4916469/files/core.tar.gz?download=1' -O db/core.tar.gz
$ cd db && tar -xvf core.tar.gz 
```

And, optionally (in the same directory):

```
$ wget 'https://zenodo.org/record/4916469/files/impute.tar.gz?download=1' -O impute.tar.gz
$ tar -xvf impute.tar.gz --strip-components=1 
```

!!! tip

    The combined size of these files is around 70GB. Make sure you have plenty
    of hard drive space in your home directory. They can take some time to
    download from zenodo.

`snpQT` is now installed. However, to run `snpQT` you will need to pick (and
possibly set up) a **profile** and make a **parameter file**.

Profiles tell `snpQT` **how** to run bioinformatics software and **where** to
run it.

The parameter file tells `snpQT` what types of things you would like to do with
your input data. It's important that you make this file and understand what each
parameter does. For convenience, we provide an example parameter file with
sensible defaults.

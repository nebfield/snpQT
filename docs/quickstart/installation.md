# Installation

Assuming you're running a modern version of Linux, firstly install
[Nextflow](https://www.nextflow.io):

```
$ curl -s https://get.nextflow.io | bash 
```

The minimal supported version of Nextflow is `21.04.3`. If you have already installed an older version of Nextflow, 
you can update the version like this:

```
$ nextflow self-update 
```

Make sure to add `nextflow` to [PATH](https://unix.stackexchange.com/a/26059) so
you can run nextflow in your terminal from anywhere. Then simply run:

```
$ git clone --branch v0.1.5 https://github.com/nebfield/snpQT.git
```

Where `--branch v0.1.5` reflects the current latest release of `snpQT`. You can find different
releases of snpQT on our [releases
page](https://github.com/nebfield/snpQT/releases).  You should see a helpful
message in your terminal if everything went well. In case you wish to pull the latest changes to the snpQT 
repository check our [Advanced Installation](https://snpqt.readthedocs.io/en/latest/quickstart/installation/) guide.

Finally, you must download some reference data
([10.5281/zenodo.4916468](https://doi.org/10.5281/zenodo.4916468)) for snpQT to
work. Download the reference files to the local path:

```
$ cd snpQT
$ mkdir db
$ wget 'https://zenodo.org/record/4916469/files/core.tar.gz?download=1' -O db/core.tar.gz
$ cd db && tar -xvf core.tar.gz 
```

And, if you wish to perform imputation you can download the reference files (in the same directory) using the following commands:

```
$ wget 'https://zenodo.org/record/4916469/files/impute.tar.gz?download=1' -O impute.tar.gz
$ tar -xvf impute.tar.gz --strip-components=1 
```

!!! tip

    The combined size of these files is around 37GB. Make sure you have plenty
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

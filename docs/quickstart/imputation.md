# Imputation set up 

The software we use for imputation, `impute5`, has a restrictive license. You
will need to accept this license, [download the
software](https://jmarchini.org/software/), and configure `snpQT` to use it.

## Singularity

To set up imputation with Singularity, firstly copy `impute5 v1.1.4` to the
`environments/impute5` directory:

```
$ cd environments/impute5/
$ cp ~/Downloads/impute5_v1.1.4.zip . 
```

Then build a Singularity `.sif` file:

```
$ sudo singularity build impute5.sif impute5.def
```

Singularity containers can run by normal users, but building Singularity images
requires **root access**. If you want to impute on a computer where you don't have
root access, you can transfer the `impute5.sif` to `environments/impute5/`.

If you don't have root access on any computers, you could try Singularity's
[remote builder service](https://cloud.sylabs.io/builder). Or maybe use a
virtual machine!

## Modules

If you are using `-profile modules` on a cluster, the easiest thing to do is to
run `module avail` and search for `impute5`. On our cluster, `impute5` is
installed as `impute/5_1.1.4`. If your cluster is different, you will need to
edit `conf/modules.config`:

```
process {
  // modules need conda too 
  conda = "$baseDir/environments/snpqt/environment.yml"  
  withName: run_snpflip {
    conda = "$baseDir/environments/snpflip/environment.yml"
  }
  withName: 'pca_prep|eigensoft' {
    conda = "bioconda::eigensoft=7.2.1 bioconda::plink=1.90b6.18"
  }
  withName: phasing {
    conda = 'bioconda::shapeit4=4.1.3'    
  }
  
  // imputation with environment modules
  // http://modules.sourceforge.net
  // https://www.nextflow.io/docs/latest/process.html#module
  // module names won't be consistent across systems :( 
  withName: 'convert_imp5|impute5' {
    module = 'impute/5_1.1.4'
  }
}
```

And change `impute/5_1.1.4` to whatever is listed in `module avail`. `-profile
modules` also requires Anaconda to be installed. 

## Docker

Docker requires root access to build and run containers, so isn't an option for
many people. However, it's simple to build a docker image from the Dockerfile:

```
$ cd environments/impute5/
$ cp ~/Downloads/impute5_v1.1.4.zip .
$ sudo docker build -t impute5 . 
```


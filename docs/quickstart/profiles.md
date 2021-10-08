# Profiles

## Where to run `snpQT`

Profiles control how `snpQT` runs. Firstly, there are two options that control
**where** modules are run:

* `-profile standard`
    * `snpQT` runs modules on a normal computer. 
* `-profile cluster`
    * `snpQT` runs modules by submitting them to a SLURM queue.

You must pick one of these profiles to run `snpQT`. Most people will use the
`standard` profile.

To configure the `cluster` profile you will need to edit the
`conf/cluster.config` file, which is in the `snpQT` local path.We have configured `cluster.config`
to run on the [Northern Ireland High Performance Centre](https://www.ni-hpc.ac.uk) by
default:

```
process {
    executor = 'slurm'
    // default queue for short jobs less than 3 hours
    queue = 'k2-hipri' 
    withLabel: bigmem {
    	memory = 128.GB
	errorStrategy = 'retry'
	maxRetries = 2
	queue = { task.attempt == 1 ? 'k2-hipri' : 'k2-medpri' } 
    }
}
```

You will need to change `k2-hipri` and `k2-medpri` to reflect the names of the
partitions in your local cluster. 
    
## How to run `snpQT`

Secondly, `snpQT` has several options to automatically install a lot of
bioinformatics software. Most people will use either `-profile conda` or
`-profile singularity`:

* `-profile conda` uses
  [Anaconda](https://docs.anaconda.com/anaconda/install/index.html) to
  automatically install software packages. Every workflow will work well except
  imputation.
* `-profile singularity` uses [Singularity](https://sylabs.io/singularity/) to
  automatically provision containers to run software packages. Every workflow
  will work well.

You will need to have installed Anaconda or Singularity for these profiles to
work.

Other profiles, including `-profile docker` and `-profile modules` are described
further in the advanced guide.

## Summary

You will need to pick a combination of profiles listed above to run `snpQT`, for
example `-profile standard,conda` or `-profile cluster,singularity`.

If you don't set a profile `snpQT` assumes that you already have installed and
configured all of the bioinformatics software we use in the pipeline, which is
quite unlikely. 

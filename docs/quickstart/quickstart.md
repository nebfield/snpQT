# Putting it all together

Assuming you've been able to:

* [Install](installation.md) `snpQT`
* Pick and maybe configure [a profile](profiles.md)
* Set up your own [parameter file](parameters.md)

You're now ready to run `snpQT`:

```
nextflow run main.nf -params-file <path_to_params_file> -profile standard,singularity -w <path_to_work_dir> -resume
```

* `<path_to_params_file>` is the filepath of a parameter file that you have
  configured e.g. `$HOME/test_run.yaml`
* `<path_to_work_dir>` is the filepath to a directory that will store
  intermediate files generated by the pipeline. This directory is used to make
  checkpoints and allows `-resume` to continue work if the pipeline stops. The
  directory can get quite big.

Once the pipeline finishes, a summary analysis and your output data will be
published to the `--results` directory you specified in the parameter file.

A simple full example using the provided toy dataset (stored in `data/` directory) including a Sample and Variant QC:

```
$ nextflow run main.nf -params-file $HOME/parameters.yaml -profile standard,conda -w $HOME/work
```

!!! Tip
	* Using a parameter file is a good practise, if you prefer using arguments on the command-line, the equivalent
	command would be: `nextflow run main.nf -profile standard,conda -params-file $HOME/parameters.yaml -w $HOME/work --qc --db db/ --bed data/toy.bed --bim data/toy.bim --fam data/toy.fam --results results/`.
 
Feel free to inspect the reports generated in `$HOME/results`. The
reports will look a bit odd because the input data is synthetic. For a more detailed exploration of all the implemented
workflows on the toy dataset, as well as a real-life ALS- Control Cohort visit our [Tutorial](https://snpqt.readthedocs.io/en/latest/user-guide/results/) page.

There are nine individual workflows in `snpQT`. We provide a detailed
explanation of each workflow in the [User guide](https://snpqt.readthedocs.io/en/latest/user-guide/background/)
section.

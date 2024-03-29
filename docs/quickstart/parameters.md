# Parameter file

The parameter file controls nearly every part of the `snpQT` pipeline, including:

* The input data location
* Which combination of workflows to run (imputation, quality control, etc.)
* Parameters for important workflow processes 
* The output data location

It's important that you understand what each parameter means. Different data
sets will require different parameters if you want to do sensible quality
control and analysis.

We recommend using a parameter file instead of specifying parameters at the
[command line](https://www.nextflow.io/docs/latest/config.html) because:

* You'll have a permanent record of your specified parameters. If you wanted to
  publish an analysis you could include this file as a supplement.
* `snpQT` has a lot of parameters. Specifying them all on the command line can
  be confusing!

To make things easier we have provided an example parameters.yaml file with
`snpQT`:

```
  # input parameters -----------------------------------------------------------
  
  # set database path to the directory of the core and impute databases
  db: 'db/'
  
  # fam filepath is mandatory for all workflows
  fam: 'data/toy.fam'

  # if you are doing build conversion, input data needs to be a VCF filepath
  # (otherwise set to false)
  vcf: false

  # if you are not doing build conversion, input data needs a .bed / .bim filepath
  bed: 'data/toy.bed'
  bim: 'data/toy.bim'
  
  # output parameters ----------------------------------------------------------
  results: 'results/'
  
  # workflow parameters --------------------------------------------------------

  # maximimum cpus per process
  cores: 4
  
  # change a workflow parameter to true, if you want to use it
  # check how you can combine multiple workflows at:
  # https://snpqt.readthedocs.io/en/latest/user-guide/background/
  
  # build conversion workflow
  convert_build: false
  # sample and variant QC workflow
  qc: true
  # population stratification workflow
  pop_strat: false 
  # local phasing and imputation workflow
  impute: false
  # pre-Imputation QC workflow
  pre_impute: false
  # post-Imputation QC workflow
  post_impute: false
  # GWAS workflow
  gwas: false
  # Download and prepare a core database passing "core"
  # or an imputation-related database passing "impute"
  # (not recommended, as it is slow, instead download the database directly from zenodo)
  download_db: false
  # help workflow
  help: false
  
  # build conversion parameters ------------------------------------------------
  
  # set to 37 or 38, if the input data are aligned in build 37 or 38,
  # accordingly
  input_build: 38
  
  # set to 37 or 38, if the output data is desired to be aligned in build 37 or
  # 38, accordingly
  output_build: 37
  
  # assign the memory size that the LiftoverVCF utility can use, in gigabytes
  mem: 16

  # qc & popstrat parameters ---------------------------------------------------
  
  # set to false if you want to skip the sex discrepancies check, 
  # recommended when your input data do not contain any sex chromosomes
  sexcheck: false
  
  # set to false if you want to remove the sex chromosomes from your dataset
  keep_sex_chroms: true
  
  # remove samples based on call rate (accepted range 0-1)
  # example below: samples with =<98% call rate are removed  
  mind: 0.02
  
  # control the pruning process
  indep_pairwise: '50 5 0.2'
  
  # remove variants based on call rate (accepted range 0-1)
  # example below: variants with =<98% call rate are removed
  variant_geno: 0.02
  
  # remove potentially related samples based on relationship-based pruning
  # threshold (accepted range 0-1)
  # example below: samples with a 3rd degree relationship and closer are removed
  king_cutoff: 0.125
  
  # change the Hardy-Weinberg Equilibrium p-value threshold (accepted range 0-1)
  hwe: 1e-7
  
  # remove variants with =< X% Minor Allele Frequency (accepted range 0-1)
  # example below: variants with =< 5% MAF are removed
  maf: 0.05
  
  # remove variants based on an X p-value threshold for missingness in
  # case/control status if you have quantitative data this step is skipped,
  # using the parameter --linear true (see below)
  missingness: 10e-7
  
  # assign population labels for the 1,000 Genome data 
  # using --popfile [super] for super population labels (e.g. EUR, AFR, AMR) 
  # or --popfile [sub] for subpopulation labels
  popfile: super
  
  # change the population codes that you wish to include in the poplist.txt file
  # that is used in smartpca
  # accepted values: --popcode [""(default), EUR/AFR/SAS... ]
  popcode: " "
  
  # change the optional parameters to the parameter file for smartpca
  # accepted values: --parfile [false (default), parfile.txt]
  parfile: false
  
  # change the default number of first Principal Components which are used to
  # create a covariates file
  # accepted range: 1-20
  pca_covars: 3
  
  # set to true if you want to remove samples with a missing phenotype
  rm_missing_pheno: false
  
  # set to false if you want to skip heterozygosity check step
  heterozygosity: true 

  # gwas parameters ------------------------------------------------------------
  
  # set to --covar_file [covar.txt], if you want to import your own covariates
  # file  
  covar_file: false
  
  # set to true if you contain quantitative data. in this case linear regression
  # will be performed  
  linear: false

  # imputation parameters ------------------------------------------------------
  
  # 128GB memory per chrom, controls the number of chromosomes that are imputed
  # at the same time.  cluster profile ignores this because it queues jobs in
  # SLURM
  # accepted range: 1-23
  impute_chroms: 1 

  # postimputation parameters --------------------------------------------------
  
  # change the info score which expresses the quality of imputation per marker
  # accepted range: 0-1
  info: 0.7
  
  # change the Minor Allele Frequency threshold in post-Imputation QC
  # accepted range: 0-1
  impute_maf: 0.01
  
  # The chosen default thresholds are used only to improve the user experience,
  # they have been chosen based on experience on our own datasets, 
  # however, each dataset is unique, so please feel free to change them
  # taking into account the accepted ranges
  
  # dummy parameters to silence nextflow warnings ------------------------------
  
  impute5_version: '_1.1.4_static'
```


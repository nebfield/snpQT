# Quickstart

(Please see the [Tutorial section](results.md) for a more detailed explanation)

[TOC]

This quickstart assumes you know what you're doing when it comes to human genomic variant data and you're not interested in [a more thorough explanation](https://snpqt.readthedocs.io/en/latest/user-guide/workflows/).

`snpQT` contains nine different workflows:

| Workflow                        | Inputs              | db             | Docker | Conda | Modules |
|---------------------------------|---------------------|----------------|--------|-------|---------|
| Reference data setup            | Not applicable      | Not applicable | Yes    | Yes   | Yes |
| Build conversion                | `--vcf --fam`       | core           | Yes    | Yes   | Yes |
| Sample quality control          | `--bed --bim --fam` | core           | Yes    | Yes   | Yes |
| Population stratification       | `--bed --bim --fam` | core           | Yes    | Yes   | Yes |
| Variant quality control         | `--bed --bim --fam` | core           | Yes    | Yes   | Yes |
| Pre-imputation quality control  | `--bed --bim --fam` | core           | Yes    | Yes   | Yes |
| Local phasing & imputation      | `--bed --bim --fam` | core + impute  | Yes    | No    | Yes |
| Post-imputation quality control | `--vcf --fam`       | core           | Yes    | Yes   | Yes |
| GWAS                            | `--bed --bim --fam` | core           | Yes    | Yes   | Yes |


The first time you run a `snpQT` workflow it may seem a little slow. Nextflow will quietly download conda environments and docker images required by the pipeline.

Below are some examples of common tasks that `snpQT` can do (after the [database set up](https://snpqt.readthedocs.io/en/latest/user-guide/installation/#core-reference-data). All examples below assume that you're currently in the `snpQT` base directory. Make sure to use different combinations of workflows depending on what is best for your needs! 

## Human genome build conversion

```nextflow
cd snpQT/
nextflow run main.nf -profile conda --vcf data/toy.vcf.gz --fam data/toy.fam --results results/ -resume --convert_build --input_build 37 --output_build 38
```

* Inputs:
    * `--vcf`: This workflow requires a valid VCF file of human genomic data
    * `--fam`: This workflow also requires an accompanying PLINK .fam file which should contain the same samples as the VCF file.
* snpQT options:
    * `--convert_build` runs the [build conversion workflow](workflows.md#build-conversion)
	* `--input_build` tells the [build conversion workflow](workflows.md#build-conversion) that the input data are aligned on b37
	* `--output_build` tells the [build conversion workflow](workflows.md#build-conversion) that the output data should be aligned on b38.
    * `--results` specifies a directory where the output files are copied to	
* Nextflow options:
    * `-resume` is helpful if you want to try different combinations of workflows later, keeping cached work
    * `-profile conda` can be replaced with `-profile docker` or `-profile modules` depending on your installation.
	
!!! note
    * `snpQT` inputs and options always start with two `--` 
    * Nextflow options always start with one `-`
    * This workflow assumes your VCF file is in human genome build 38 and it will be converted to b37. This means that in this case you can avoid using the parameters `--input_build 38 --output_build 37`
	* You can convert in the end of the `snpQT` pipeline your dataset back to b38 using `--input_build 37 --output_build 38`.

!!!Warning
	The .fam file should contain the same samples as the VCF file. It is used to update the missing sex and phenotype status of the VCF file.
	
## Quality control

```nextflow
cd snpQT/
nextflow run main.nf -profile conda --bed data/toy.bed --fam data/toy.fam --bim data/toy.bim --results results/ --qc --sexcheck false -resume
```

* Inputs:
    * `--bed`, `--fam` and `--bim` expect valid PLINK .bed, .bim and .fam files, respectively. Learn more [here](https://www.cog-genomics.org/plink2/formats).
* snpQT options:
    * `--qc` runs the [quality control workflow](https://snpqt.readthedocs.io/en/latest/user-guide/workflows/#quality-control)
    * `--results` specifies a directory where output files are copied to
    * `--sexcheck false`, the toy dataset doesn't have sex chromosomes so to avoid PLINK producing an error, `snpQT` skips the sex discrepancies check step using this parameter. 
* Nextflow options:
    * `-resume` is helpful if you want to try different combinations of workflows later
    * `-profile conda` can be replaced with `-profile docker` or `-profile modules` depending on your installation.

!!! note
    * This workflow assumes your genomic data are in human genome build 37 
	* You can combine this workflow with `--convert_built` workflow. In this case, the input files should be a VCF and a .fam file.
	
## Population stratification with quality control

```nextflow
nextflow run main.nf -profile conda --bed data/toy.bed --bim data/toy.bim --fam data/toy.fam --qc --pop_strat --sexcheck false -resume --results results_popStrat/
```

* Inputs:
    * `--bed`, `--fam` and `--bim` expect valid PLINK .bed, .bim and .fam files, respectively. Learn more [here](https://www.cog-genomics.org/plink2/formats).
* snpQT options:
    * `--qc` runs the [quality control workflow](https://snpqt.readthedocs.io/en/latest/user-guide/workflows/#quality-control)
    * `--pop_strat`runs the [population stratification workflow](workflows.md#population-stratification)
    * `--sexcheck false`, the toy dataset doesn't have sex chromosomes so to avoid PLINK producing an error, `snpQT` skips the sex discrepancies check step using this parameter
	* `--results` specifies a directory where output files are copied to.
	
!!! note
    * This workflow assumes your binary PLINK files are already aligned in b37 
	* The SNP IDs need to be in a "rs" ID format, in order to be merged with the 1,000 human genome data
	* Population stratification can not be used without the quality control workflow. If you forget to add `--qc`, `snpQT` will produce an error and will remind you to add the aforementioned workflow.
	* You can combine this workflow with `--convert_built` workflow. In this case, the input files should be a VCF and a .fam file.
	

## GWAS

```nextflow
nextflow run main.nf -profile conda --bed data/toy.bed --bim data/toy.bim --fam data/toy.fam --qc --pop_strat --gwas --sexcheck false -resume
```

* Inputs:
    * `--bed`, `--fam` and `--bim` expect valid PLINK .bed, .bim and .fam files, respectively. Learn more [here](https://www.cog-genomics.org/plink2/formats).
* snpQT options:
    * `--qc` runs the [quality control workflow](https://snpqt.readthedocs.io/en/latest/user-guide/workflows/#quality-control)
    * `--pop_strat`runs the [population stratification workflow](workflows.md#population-stratification)
    * `--gwas` runs the [GWAS workflow](workflows.md#genome-wide-association-study)
    * `--sexcheck false`, the toy dataset doesn't have sex chromosomes so to avoid PLINK producing an error, `snpQT` skips the sex discrepancies check step using this parameter.

!!! note
    * This workflow requires the `--qc` workflow in order to run. But, do not worry, if you forget to add `--qc`, `snpQT` will produce an error and will remind you to add the aforementioned workflow.
	* If your input data are built on b38 you can add `--convert_build --input_build 38 --output_build 37`. In this case, the input files should be a VCF and a .fam file.
	
## Pre-imputation quality control

```nextflow
nextflow run main.nf -profile conda --bed data/toy.bed --bim data/toy.bim --fam data/toy.fam --qc --pop_strat --pre_impute --sexcheck false -resume
```

* Inputs:
    * `--bed`, `--fam` and `--bim` expect valid PLINK .bed, .bim and .fam files, respectively. Learn more [here](https://www.cog-genomics.org/plink2/formats).
* snpQT options:
    * `--qc` runs the [quality control workflow](https://snpqt.readthedocs.io/en/latest/user-guide/workflows/#quality-control)
    * `--pop_strat`runs the [population stratification workflow](workflows.md#population-stratification)
    * `--pre_impute` runs the [pre-imputation workflow](https://snpqt.readthedocs.io/en/latest/user-guide/workflows/#pre-imputation-quality-control) preparing your dataset for phasing and imputation
    * `--sexcheck false`, the toy dataset doesn't have sex chromosomes so to avoid PLINK producing an error, `snpQT` skips the sex discrepancies check step using this parameter.

!!! note
    * This workflow assumes your VCF file is aligned in human genome build 37
	* `--pre_impute` can not be combined with GWAS workflow. It is designed for users who wish to upload a prepared clean genomic dataset to an external imputation server, or to perform imputation using a reference panel other than the latest release of the 1,000 genome data that `snpQT` supports
	* Pre-imputation workflow is nested under the `--impute` workflow (below), so `--pre_impute` and `--impute` should not be used together (when used an error message is displayed to the standard output).

## Imputation 
    
```nextflow
nextflow run main.nf -profile docker --bed data/toy.bed --bim data/toy.bim --fam data/toy.fam --qc --pop_strat --impute --sexcheck false -resume
```

* Inputs:
    * `--bed`, `--fam` and `--bim` expect valid PLINK .bed, .bim and .fam files, respectively. Learn more [here](https://www.cog-genomics.org/plink2/formats).
* snpQT options:
    * `--qc` runs the [quality control workflow](https://snpqt.readthedocs.io/en/latest/user-guide/workflows/#quality-control)
    * `--pop_strat` runs the [population stratification workflow](workflows.md#population-stratification)
    * `--impute` runs the [pre-imputation](https://snpqt.readthedocs.io/en/latest/user-guide/workflows/#pre-imputation-quality-control), [imputation](workflows.md#phasing-imputation) and [post-imputation quality control](https://snpqt.readthedocs.io/en/latest/user-guide/workflows/#post-imputation-quality-control) workflows
    * `--sexcheck false`, the toy dataset doesn't have sex chromosomes so to avoid PLINK producing an error, `snpQT` skips the sex discrepancies check step using this parameter.

!!! note
	* You will need to have set up a [docker imputation image](https://snpqt.readthedocs.io/en/latest/user-guide/installation/#imputation-prerequisite) and [additional reference database](installation.md#imputation-reference-data)
    * This workflow assumes your VCF file is aligned to human genome build 37. If your input data are built on b38 you can add `--convert_build --input_build 38 --output_build 37`. In this case, the input files should be a VCF and a .fam file.
	* Imputation workflow combines pre-imputation and post-imputation, so `--pre_impute`, `--post_impute` and `--impute` should not be used together (if they do an error message is displayed to the standard output).



## Imputation with GWAS

```nextflow
nextflow run main.nf -profile docker --bed data/toy.bed --bim data/toy.bim --fam data/toy.fam --qc --pop_strat --impute --gwas --sexcheck false -resume
```
* Inputs:
    * `--bed`, `--fam` and `--bim` expect valid PLINK .bed, .bim and .fam files, respectively. Learn more [here](https://www.cog-genomics.org/plink2/formats).
* snpQT options:
    * `--qc` runs the [quality control workflow](https://snpqt.readthedocs.io/en/latest/user-guide/workflows/#quality-control)
    * `--pop_strat`runs the [population stratification workflow](workflows.md#population-stratification)
    * `--gwas` runs the [GWAS workflow](workflows.md#genome-wide-association-study)
    * `--impute` runs the [pre-imputation](https://snpqt.readthedocs.io/en/latest/user-guide/workflows/#pre-imputation-quality-control), [imputation](workflows.md#phasing-imputation) and [post-imputation quality control](https://snpqt.readthedocs.io/en/latest/user-guide/workflows/#post-imputation-quality-control) workflows
    * `--sexcheck false`, the toy dataset doesn't have sex chromosomes so to avoid PLINK producing an error, `snpQT` skips the sex discrepancies check step using this parameter.

!!! note
    This workflow assumes your VCF file is aligned to human genome build 37. If your input data are aligned in b38 you can add `--convert_build --input_build 38 --output_build 37`. In this case, the input files should be a VCF and a .fam file.

## Post-imputation quality control

```nextflow
cd snpQT/
nextflow run main.nf -profile conda --vcf data/toy.vcf.gz --fam data/toy.fam --results results/ -resume --post_impute --sexcheck false
```

* Inputs:
    * `--vcf`: This workflow requires a valid VCF file of human genomic data
    * `--fam`: This workflow also requires an accompanying PLINK .fam file which should contain the same samples as the VCF file.
* snpQT options:
	* `--post_impute` runs the [post-imputation quality control workflow](https://snpqt.readthedocs.io/en/latest/user-guide/workflows/#post-imputation-quality-control)
	* `--sexcheck false`, the toy dataset doesn't have sex chromosomes so to avoid PLINK producing an error, `snpQT` skips the sex discrepancies check step using this parameter.

!!! note
    * `--post_impute` is designed for users who have imputed their data and they wish to perform a post-imputation quality control. You can use the output binary files of this workflow to run any other workflow that is supported by `snpQT`.
	* Post-imputation workflow is nested under the `--impute` workflow, so `--post_impute` and `--impute` should not be used together (when used an error message is displayed to the standard output).
	
## Help page

```nextflow
cd snpQT/
nextflow run main.nf --help
```

For a quick visit in the basic workflows and optional parameters, navigate through the help page in the command line.

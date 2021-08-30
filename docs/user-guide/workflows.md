# Workflow descriptions

[TOC]

## Build conversion

```nextflow
--convert_build

Build conversion workflow options:
  Mandatory:
    --vcf                                   
    --fam   
	
  Optional:	
	--input_build [38 (default),37]
	--output_build [37 (default),38]
	--mem [16 (default)]
```

`snpQT` assumes that your genomic data are built on human genome build 37. Even though b37 is not the most recent build, it is the most frequently used build among current public, reference genomic datasets (e.g. 1,000 Genome data, Haplotype Reference Concortium panel) and available SNP-array datasets, and for this reason we decided that `snpQT` should support b37 for the most of the workflows. However, if you wish to impute your data using an external online server that uses a reference panel in b38 like TOPMed, we provide input parameters to convert your `snpQT` clean dataset from b37 back to b38 (using `--input_build 37` and `--output_build 38`). 

This workflow is intented for those users whose data are built on b38 (default) or b37 and they wish to convert them to b37 or b38, accordingly, using `Picard`'s `LiftoverVcf` utility, so that their data are compatible with subsequent `snpQT` workflows and also the users can convert their data in their initial build if they wish to. If your data are built neither on b38 nor b37 then visit our [public git repo](https://github.com/ChristinaVasil/Quality-Control-Pipeline-for-Genomic-data/tree/master/1_HumanGenomeBuildConversion) where you can run the code of this workflow yourselves, using alternative .chain.gz files depending on your current build. Lastly, `--mem` parameter controls the memory size that the `LiftoverVCF` utility can use.


!!! Warning
	Be aware that a small portion of SNPs is expected to be removed between conversion. A portion of SNPs are different among different hgb versions (e.g. not present, merged or relabelled).


## Quality control

```nextflow
--qc 

Quality control workflow options:
  Mandatory:
      --bed                                   
      --bim                                   
      --fam                                   
  Optional:
    --mind  [0.02 (default), 0-1]
    --indep_pairwise ["50 5 0.2" (default), ""]
    --variant_geno [0.02 (default), 0-1]
    --hwe  [1e-7 (default), 0-1]
    --maf [0.05 (default), 0-1]
    --missingness [1e-7(default), 0-1]
	--sexcheck [true (default),false]
	--rm_missing_pheno [false (default),true]
	--keep_sex_chroms [true (default),false]
	--heterozygosity [true (default),false]
	--king_cutoff [0.125 (default), 0-1]
	--pca_covars [3 (default), 1-20]
	--linear [false (default),true]
```

### Sample quality control

Below we list the checks which are followed in the Sample QC workflow:

- **Missing variant call rate check**: Remove very poor quality SNPs based on call rate (`snpqt` default threshold is 0.1). This way we avoid removing samples based on very poor quality SNPs (which we would remove later in Variant QC workflow anyway). 

- **Missing sample call rate check**: Remove samples with lower than 98% call rate (default threshold `--mind 0.02`) and then visualize the distribution for all samples using histograms and scatterplots before and after the applied threshold. This threshold can be changed using the `--mind` parameter.

- **Check for sex discrepancies**: Remove problematic samples for which (1) pedigree sex does not agree with the predicted sex based on sex chromosome homozygosity or (2) there is no sex phenotype.  This step is important to avoid potential sample mix-ups and/or DNA contaminations.

!!! Warning
	This step will not run if sex chromosomes are not included in the .bim file. In this case, use the `--sexcheck false` parameter to skip this step.

- **Removal of non-autosomal SNPs**: The default mode of `snpQT` is to keep the sex chromosomes. If the user wishes to remove the sex chromosomes, use `--keep_sex_chroms false` 

- **Heterozygosity check**: Identify and remove heterozygosity outliers (samples that deviate more than 3 units of Standard Deviation from the mean heterozygosity). The distribution of the samples' heterozygosity is vizualized through a histogram and scatterplot. Extreme heterozygosity implies inbreeding and/or DNA contamination. This step can be skipped using `--heterozygosity false`. This can be useful in occasions where the dataset has been already through heterozygosity once, so the user does not wish to remove samples again based on this metric.

- **Check for cryptic relatedness and duplicates**: Check for cryptic pairs of relatives using `plink2`'s [relationship-based pruning](https://www.cog-genomics.org/plink/2.0/distance#king_cutoff) threshold. Relatedness is defined as 3rd degree or closer (default threshold for this step is 0.125). This threshold can be changed using the `--king_cutoff` parameter.

- **Removal of samples with a missing phenotype**: Remove samples with missing phenotypes. As *missing phenotype* here we refer to case/control status (i.e. the last column in your PLINK .fam file). The default option in `snpQT` is to skip this step.

At the end of Sample QC a .log file is generated listing the number of samples, variants, phenotypes and working directory for each step where the intermediate files are stored, as well as an .html report containing before-and-after the chosen thesholds plots for the vast majority of the aforementioned steps.

### Variant quality control

The Variant QC workflow is the second part of the `–-qc` utility of `snpQT`. It is good practice to first filter low quality samples in order to reduce the risk of removing a potentially high-risk variant during Variant QC. For this reason, the population stratification workflow (if chosen to run by the user), which is essentially a Sample QC step, is designed to run in between Sample QC and Variant QC. Below we list the checks which are followed in the Variant QC workflow:

- **Missing variant call rate check **: Remove poor quality SNPs based on a more strict call rate threshold (`snpqt` default threshold is 0.02). This threshold can be changed using the parameter `–-variant_geno`.

- **Hardy-Weinberg equilibrium (HWE) deviation check**: Remove SNPs that significantly deviate from the Hardy-Weinberg equilibrium (HWE) (`snpqt` default threshold is hwe p-value < 10e-7), indicating a genotyping error, and visualize the distribution of SNPs with extreme deviation.  This threshold can be changed using the parameter `–-hwe`.

- **Minor Allele Frequency (MAF) check**: Remove SNPs with low MAF (`snpqt` default threshold is 0.05) and visualize the MAF distribution. This threshold can be changed using the `–-maf` paramater in `snpQT`. Rare SNPs (having a very low MAF) are usually considered as false-positives and need to be excluded from further analysis.

- **Missingness in case/ control status check**: Remove SNPs with a statistically significant association of missingness (low call rate) and case/control status (`snpqt` default threshold is p-value < 10e-7). The threshold can be changed using the parameter `–-missingness`.

!!!Warning
	Missingness in case/ control status check can not be performed in quantitative data. If you do not have binary data use the `--linear` parameter to skip this step.
	
- **Generate covariates using the first X Principal Components of each sample**: `snpQT` by default uses the first 3 Principal Components (PCs) to account for inner population structure. The number of PCs can by changed using the `--pca_covars` parameter which can take as input a number from 1 to 20, with 1 starting from the first Principal Component of the PCA. Output from this step is used in the `--gwas` workflow, in order to account for potential inner population structure (it is best used along with `--pop_strat`).

At the end of Variant QC we provide an .html report containing the following: 

- Before-and-after applying the chosen thesholds plots for all the aforementioned steps

- Three 2D PCA plots (PC1vsPC2, PC1vsPC3 and PC2vsPC3) of the user's data annotated with case/control status

- 3D interactive PCA plot of only the user's data annotated with case/control status.

We also provide all the figures, logs and binary plink files in separate directories `./qc/figures/`, `./qc/logs/` and `./qc/bfiles/`, respectively. 

## Population stratification

```nextflow
--pop_strat

Population stratification options:
  Mandatory:
    --qc
  
  Optional:
    --variant_geno [0.02 (default), 0-1]
    --indep-pairwise ["50 5 0.2" (default), ""]
    --popfile [super (default), sub]
	--parfile [false (default), parfile.txt]
	--popcode [""(default), "EUR"/"AFR"/"SAS"... ]
```

This workflow aims to identify and remove outliers using [EBI's phased latest release 1,000 human genome reference panel](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) aligned in human genome build 37, as well as account for potential inner population structure in the later GWAS workflow.  Population stratification is an essential step in QC analysis, since it minimizes the possibility that the difference in the allele frequencies is caused by the different ancestry of the samples.

!!!Note
	Population stratification is a sample QC step and for this reason it is designed to run between the Sample and Variant QC `--qc` workflows.

The first step in `--pop_strat` is to prepare and merge 1,000 human genome data with the user's dataset following a number of QC steps using `plink2`. The processing steps for 1,000 human genome data include:

- Removing sex chromosomes

- Checking that all alleles are on the forward strand

- Removing duplicated SNPs

- Removing multi-allelic variants

- Removing very poorly genotyped variants

- Removing samples with low call rate (>98%)

- Removing poorly genotyped variants in a second more stringent threshold (`--variant_geno 0.02`)

- Removing rare variants

- Keeping only highly independent SNPs by excluding regions with high linkage equilibrium and pruning using `--indep-pairwise 50 5 0.2`.

The next step is to prepare the user's dataset (samples having already been processed through the Sample QC workflow):

- Removing poorly genotyped variants in a second more stringent threshold (`--variant_geno 0.02`)

- Removing rare variants

- Keeping only highly independent SNPs by excluding regions with high linkage equilibrium and pruning using `--indep-pairwise 50 5 0.2`.

- Checking that all alleles are on the forward strand

- Removing ambiguous SNPs

- Removing duplicated SNPs

When both datasets are prepared, the next step is merging: keeping only mutual SNPs shared by both the user's dataset and the 1,000 human genome data. 

Then we make a popfile summarizing the samples of the merged dataset adding a third column labelling the population origin of each sample. User's samples are automatically labelled as "OWN". The population label for 1,000 human genome data is defined by the `--popfile` flag, having as a default to use super population labels (e.g. EUR, AFR, AMR). If you want to use subpopulation labels then you should add the `--popfile sub` parameter.

When the popfile and the merged dataset are ready, it is time to run Eigensoft's `smartpca` software. `smartpca` needs a set of parameters in order to run, which are in the form of a file (parfile). We provide the option to change this parfile according to the users' needs (you can have a look at the parameters of the parfile [here](https://github.com/DReichLab/EIG/tree/master/POPGEN)). Our default parameters for the parfile are the following:

```
	genotypename: C6_indep.bed
	snpname:      C6_indep.bim
	indivname:    C6_indep.pedind
	evecoutname:  eigenvec
	evaloutname:  eigenval
	outlieroutname: excluded_outliers.txt
	numthreads:   2
	poplistname: poplist.txt
	numoutlierevec: 7
	autoshrink: YES
```

The user can change/add/remove parameters by passing a file to the `--parfile [file]` parameter. However, the following parameters are mandatory for `snpQT` to run, and therefore should not be altered or removed by the user:

```
	genotypename: C6_indep.bed
	snpname:      C6_indep.bim
	indivname:    C6_indep.pedind
	evecoutname:  eigenvec
	evaloutname:  eigenval
	outlieroutname: excluded_outliers.txt
```

When `--pop_strat` has finished, we provide the following PCA results in a .html report:

- 3D interactive PCA plot of merged dataset containing all samples **before and after outlier removal**

- 2D PCA plots of merged dataset containing all samples

The last and most important step in `--pop_strat` is the:

- **Removal of any samples that are outliers**. 

We also provide all the figures, logs and binary plink files in separate directories `./pop_strat/figures/`, `./pop_strat/logs/` and `./pop_strat/bfiles/`, respectively. 

## Pre-Imputation quality control

```nextflow
--pre_impute

Pre-imputation workflow options:
    Mandatory:
      --qc

```


In pre-imputation workflow the user's dataset is prepared for phasing and imputation, either locally using `snpqt` or for the purpose to upload a clean vcf.gz file to an external imputation server. 

The main aims of this workflow are:
		
- **Flip SNPs that are on the reverse strand**
- **Remove ambiguous SNPs**
- **Remove one of each pair of duplicated SNPs**
- **Fix and swap the reference allele** using `bcftools` `+fixref` plug-in.

A clean and properly prepared for imputation vcf file and a index file are stored in `./pre_imputation/files/` directory.

!!!Warning
	* Pre-Imputation workflow can not be combined with `--gwas`, as its purpose is to upload a properly prepared vcf file to an external imputation server
	* Pre-imputation and post-imputation workflows are run internally when `--impute` is used to run imputation locally. So, you do not have to combine `--impute` with `--pre_impute` and `--post_impute` workflow parameters. If you do, `snpQT` is designed to throw an error.

## Phasing & Imputation

```nextflow
--impute

Imputation workflow options:
    Mandatory:
      --qc
	 
	 Optional:
	  --impute_maf [0.01 (default), 0-1]
	  --info [0.7 (default), 0-1]
	  --impute_chroms [1(default), 1-23]

```

This workflow has three main parts:

- **Phasing**: Perform phasing using `shapeit4` and index all phased chromosomes

- **Imputation**: Perfom local imputation using `impute5`, following the steps below:
	
	* Convert reference genome (1,000 human genome) into a `.imp5` format for each chromosome, using 	imp5Converter
	* Run `impute5` using the converted reference genome, genetic maps and user's prepared phased data.
	* Using the parameter `--impute_chroms [1(default), 1-23]` you can control the number of chromosomes that are imputed at the same time. As higher the number of chromosomes is, the more RAM your machine will need to use.
		
## Post-imputation quality control

```nextflow
--post_impute

Post-imputation workflow options:
    Mandatory:
	  --vcf
	  --fam
	
	Optional:
	  --impute_maf [0.01 (default), 0-1]
	  --info [0.7 (default), 0-1]

```
Tha main aims of this workflow are:

- **Merge all imputed chromosomes with `bcftools`**

!!! Warning
	We use `-n` parameter during concatenation which makes the process run much faster but assumes that your .vcf.gz files are sorted.
	
- **Filter variants**: We filter all poorly imputed variants based on info score (default is 0.7 - if you want to change the threshold use `--info`) and filter based on MAF (default is 0.01 - if you want to change the threshold use `--maf`)

- **Annotate missing SNP ids**: The annotation of the missing SNPs is in this format Chromosome:Position:Reference_Allele:Alternative_Allele

- **Handle all categories of "duplicated" SNP ids**: 

	- Remove exact duplicates
	
	- Remove multi-allelics
	
	- Annotate merged SNPs: Merged SNPs are a special category of duplicated SNP ids with different position and ref/alt alleles. You can read more about this category [here](https://www.ncbi.nlm.nih.gov/books/NBK44468/#!po=6.25000).

!!! Note
	Even though Post-Imputation workflow is nested under the `--impute` flag, it is also designed to run independently as some users might prefer running imputation using online servers, or have already imputed data and they wish to proceed with a Post-Imputation QC. 
	
!!! Warning
	The input `--vcf input.vcf.gz` file should contain the same samples as the input `--fam input.fam` file. In case, they contain different samples `snpQT` will output an error.
	
## Genome-Wide Association Study

```nextflow
--gwas

GWAS workflow options:
    Mandatory:
      --qc

    Optional:
	  --pop_strat
      --impute  
	  --covar_file [false (default), covar.txt]
	  --pca_covars [3 (default), 1-20]
	  --linear [false (default),true]
```

The Genome-Wide Association Studies (GWAS) workflow aims to identify markers with a statistically significant association with the trait of interest. This workflow performs logistic or linear regression (depending on the phenotype) with and without covariates (calculated at the end of the `---qc` workflow and passed to the GWAS workflow). The `---gwas` workflow illustrates the results of Generalized Linear Regression (`plink2`'s `--glm`) in the forms of a Manhattan plot and a Q-Q plot. The main processes of this workflow include:

- **Run logistic or linear regression**: 

	* Adjusting for covariates accounting for fine-scale population structure. The covariates can either be imported by the user using the `--covar_file covar.txt` parameter or can be generated in `--qc` workflow (for better results combine with `--pop_strat`) using the first X Principal Components of the generated PCA using the `--pca_covars [3 (default), 1-20]` parameter.
	
	* Not adjusting for covariates.
	
!!!Warning 
	* The `--covar_file` parameter can not be combined with `--pca_covars`.
	* For the format of the covar.txt file please advise [plink2](https://www.cog-genomics.org/plink2/).
	* Since we are using `plink2`'s `--glm`, the covariates file can not contain columns for the sex of the samples. Sex is automatically accounted for in `--glm`.
	
We added these two processes for two main reasons. The first one is to make `--gwas` useful for users who do not wish to run `--pop_strat` or use covariates or even inspect/compare the effect of the used covariates. In this case, the first process will not produce an output (designed so as to not produce an error), but the second process will provide the expected results (along with a Manhattan plot and a Q-Q plot). The second reason is for users that wish to run `--pop_strat`, use `--pop_strat` covariates or insert their own covariates and it would be helpful for them to compare their GWAS results in the output plots with and without covariates.

- **Illustrate a Q-Q (Quantile-Quantile) plot**: A plot to inspect the lambda genomic inflation parameter (calculated as median 1df chi-square stat / 0.456), expressing the relationship between the observed and the expected quantile p-values of the sample cohort under a normal distribution.

- **Illustrate Manhattan plot**: A plot where the association p-values are represented on the y-axis, and the positions of the tested variants are represented on the x-axis.

## Log steps table 

The table below links every step with a log code which can be useful for users who want to inspect the log plots. The logs are generated automatically and show the removed samples and variants for each step.

| Code     | Step                                       | Flag |  Workflow|
|----------|--------------------------------------------|------|----------|
| A1       | Download and process reference files |`--download_db`| Database set up|
| A2       | Decompress reference files |`--download_db`| Database set up|
| A3       | Remove duplicate records in reference files |`--download_db`| Database set up|
| A4       | Download and unzip genetic maps for phasing |`--download_db`| Database set up|
| A5       | Annotate variants to a "chr:pos:ref:alt" format |`--download_db`| Database set up|
| B1       | Decompress FASTA file                             | `--convert_build` | Build Conversion |
| B2       | Create a dictionary file                          | `--convert_build` | Build Conversion |
| B3       | Change the chromosome ids                         | `--convert_build` | Build Conversion |
| B4       | Run LiftOver to map genome build                  | `--convert_build` | Build Conversion |
| B5       | Change the chromosome ids                         | `--convert_build` | Build Conversion |
| B6       | Convert VCF to PLINK format and update phenotypes | `--convert_build` | Build Conversion |
| C1       | Missing variant call rate check            |`--qc`|Sample QC|
| C2       | Missing sample call rate check             |`--qc`|Sample QC |
| C3       | Check for sex discrepancies                |`--qc`|Sample QC |
| C4       | Removal of non-autosomal SNPs              |`--qc`|Sample QC |
| C5       | Heterozygosity check                       |`--qc`|Sample QC |
| C6       | Check for cryptic relatedness              |`--qc`|Sample QC |
| C7       | Removal of samples with missing phenotype  |`--qc`|Sample QC |
| E8       | Missing variant call rate check            |`--qc`|Variant QC |
| E9       | Hardy-Weinberg equilibrium (HWE) deviation check|`--qc`|Variant QC |
| E10      | Minor Allele Frequency (MAF) check         |`--qc`|Variant QC |
| E11      | Missingness in case/ control status check  |`--qc`|Variant QC |
| E12      | Create a covariate file for GWAS and perform PCA |`--qc`|Variant QC|
| E13      | Combine all .log files into one file  |`--qc`| Variant QC, Sample QC |
| E14      | Create a .html report  |`--qc`| Variant QC, Sample QC |
| D3       | QC and preparation of user's data |`--pop_strat`| PopStrat|
| D4       | Fix strand errors and remove ambiguous SNPs |`--pop_strat`| PopStrat|
| D5       | Align reference allele according to reference genome |`--pop_strat`| PopStrat|
| D6       | Merge user's dataset with 1,000 Human Genome data |`--pop_strat`| PopStrat|
| D7       | Create a population file |`--pop_strat`| PopStrat|
| D8       | Principal Component Analysis & `smartpca` |`--pop_strat`| PopStrat|
| D9       | Remove outliers from user's dataset |`--pop_strat`| PopStrat|
| F1       | Set appropriate chromosome codes |`--impute`, `--pre_impute`| Pre-Imputation|
| F2/D4    | Check for strand issues and remove ambiguous SNPs|`--impute`, `--pre_impute`| Pre-Imputation|
| F3       | Remove one of each pair of duplicated SNPs  |`--impute`, `--pre_impute`| Pre-Imputation|
| F4       | Convert Plink file into BCF |`--impute`, `--pre_impute`| Pre-Imputation|
| F5       | Check and fix the REF allele `bcftools +fixref` |`--impute`, `--pre_impute`| Pre-Imputation|
| F6       | Sort the BCF, Convert `.bcf` file to `.vcf.gz` and index |`--impute`, `--pre_impute`| Pre-Imputation|
| A6       | Annotate user's variants to a "chr:pos:ref:alt" format | `--impute`| Imputation|
| G1       | Split vcf.gz file in chromosomes and index them |`--impute`| Imputation|
| G2       | Perform phasing using `shapeit4`  |`--impute`| Imputation|
| G3       | Index phased chromosomes |`--impute`| Imputation|
| G4       | Tabix reference files |`--impute`| Imputation|
| G5      | Convert vcf.gz reference genome chromosomes into .imp5 format|`--impute`| Imputation|
| G6       | Perform imputation using impute5 |`--impute`| Imputation|
| G7       | Merge all imputed chromosomes |`--impute`| Imputation|
| E13      | Combine all .log files into one file  |`--impute`| Imputation |
| H1       | Filter all poorly imputed variants using info score |`--impute`, `--post_impute`| Post-Imputation|
| H2       | Filter imputed variants using MAF |`--impute`, `--post_impute`| Post-Imputation|
| H3       | Remove duplicate variants |`--impute`, `--post_impute`| Post-Imputation|
| H4       | Remove multi-allelics variants |`--impute`, `--post_impute`| Post-Imputation|
| H5       | Identify merged variants |`--impute`, `--post_impute`| Post-Imputation|
| H6       | Update sample ids information |`--impute`, `--post_impute`| Post-Imputation|
| H7       | Update phenotype information |`--impute`, `--post_impute`| Post-Imputation|
| E13      | Combine all .log files into one file  |`--impute`, `--post_impute`| Imputation |
| I1       | Perform logistic regression with and without covariates and plot results|`--gwas`| GWAS|
| E13      | Combine all .log files into one file  |`--gwas`| GWAS |
| E14      | Create a .html report  |`--gwas`| GWAS |

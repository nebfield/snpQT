def printHelp() {
  log.info"""
  Usage:
    nextflow run main.nf -profile (docker,conda) ( --download_db | --convert_build | --qc | --pop_strat | --impute | --pre_impute | --post_impute | --gwas ) [workflow-options]
  Description:
    Make your single-nucleotide polymorphisms cute with snpQT. A pipeline for quality control, population stratification, imputation, and GWAS of human genomic variants.
  
  Full documentation available at: https://snpqt.readthedocs.org

  Nextflow arguments:
    -profile 
	-resume
  
  Main workflow arguments:
    --download_db                                           Run the reference database setup workflow
    --convert_build                                         Run the build conversion workflow
    --qc                                                    Run the quality control workflow
    --pop_strat                                             Run the population stratification workflow
    --pre_impute                                            Run the pre-imputation quality control workflow
	--impute                                                Run the imputation workflow
    --post_impute                                           Run the post-imputation quality control workflow
    --gwas                                                  Run the genome wide association study workflow
   
    --results                                               The output directory for results files (default: snpQT/results)

  Database download workflow options:
    --download_db core                                      Download the core reference files to enable --convert_build, --qc, --pop_strat, and --gwas
    --download_db impute                                    Download additional reference files to enable --impute

  Build conversion workflow options:
    Mandatory:
		--vcf                                               Path to VCF file (b38 or b37)
		--fam                                               Path to fam file containing phenotypes

	Optional:	
		--input_build [38 (default),37]                     The human genome build that input vcf is aligned
		--output_build [37 (default),38]                    The human genome build that output vcf is aligned
	
  Quality control workflow options:
    Mandatory:
		--bed                                               Path to bed file
		--bim                                               Path to bim file
		--fam                                               Path to fam file
   
    Optional:
		--mind  [0.02 (default), 0-1]                       Threshold for missing sample call rate 
		--indep_pairwise ["50 5 0.2" (default), ""]         Parameters for identifying independent markers
		--variant_geno [0.02 (default), 0-1]                Threshold for missing variant call rate
		--hwe  [1e-7 (default), 0-1]                        Threshold for Hardy-Weinberg Equilibrium deviation test
		--maf [0.05 (default), 0-1]                         Threshold for minor allele frequency 
		--missingness [1e-7(default), 0-1]                  Threshold for the association between missingness and phenotype status
		--sexcheck [true (default),false]                   If false, skip the check for sex discrepanies (use if your input genomic data do not contain sex chromosomes)
		--heterozygosity [true (default),false]             If false, skip the heterozygosity pruning (still able to inspect the distribution)
		--keep_sex_chroms [true (default),false]            If false, sex chromosomes are removed
		--pihat [0.125 (default), 0-1]                      Threshold for relatedness identification
		--rm_missing_pheno [false (default),true]           Remove samples with a missing phenotype
		--pca_covars [3 (default), 1-20]                    Number of Principal Components to be accounted for generating covariates
	    --linear [false (default),true]                     Use if you have quantitative data
    
  Population stratification options:
    Mandatory:
		--qc                                               Run the quality control workflow

    Optional:
		--indep-pairwise ["50 5 0.2" (default), ""]        Parameters for identifying independent markers
		--racefile [super (default), sub]                  If set to super or sub, superpopulation and subpopulation codes are used
		--parfile [false (default), parfile.txt]           Add a parfile for smartpca
		--racecode [""(default), "EUR"/"AFR"/"SAS"... ]    List population codes for smartpca

  Pre-imputation workflow options:
    Mandatory:
		--qc                                               Run the quality control workflow
		--pre_impute                                       Run the pre-imputation quality control workflow
	
  Imputation workflow options:
    Mandatory:
		--qc                                               Run the quality control workflow

    Optional:
		--impute_maf [0.01 (default), 0-1]                 Threshold for minor allele frequency 
		--info [0.7 (default), 0-1]                        Threshold for info score
	
  Post-imputation workflow options:
    Mandatory:
		--vcf                                              Path to VCF file (b38 or b37)
		--fam                                              Path to fam file containing phenotypes
	
    Optional:
		--impute_maf [0.01 (default), 0-1]                 Threshold for minor allele frequency 
		--info [0.7 (default), 0-1]                        Threshold for info score

  GWAS workflow options:
    Mandatory:
		--qc                                               Run the quality control workflow

    Optional:
		--pop_strat                                        Run the population stratification workflow
		--impute                                           Run the imputation workflow
		--covar_file [false (default), covar.txt]          Import a covariates file (can not be used along with --pca_covars)
		--pca_covars [3 (default), 1-20]                   Number of Principal Components to be accounted for generating covariates
		--linear [false (default),true]                    Use if you have quantitative data 
  """.stripIndent()
}
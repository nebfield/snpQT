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
    --download_db                             Run the reference database setup workflow
    --convert_build                           Run the build conversion workflow
    --qc                                      Run the quality control workflow
    --pop_strat                               Run the population stratification workflow
    --pre_impute                              Run the pre-imputation quality control workflow
	--impute                                  Run the imputation workflow
    --post_impute                             Run the post-imputation quality control workflow
    --gwas                                    Run the genome wide association study workflow
   
    --results                                 The output directory for results files (default: snpQT/results)

  Database download workflow options:
      --download_db core                      Download the core reference files to enable --convertBuild, --qc, --popStrat, and --gwas
      --download_db impute                    Download additional reference files to enable --impute

  Build conversion workflow options:
    Mandatory:
      --vcf                                   Path to VCF file (b38 or b37)
      --fam                                   Path to fam file containing phenotypes

	Optional:	
	  --input_build [38 (default),37]
	  --output_build [37 (default),38]
	
  Quality control workflow options:
    Mandatory:
      --bed                                   Path to bed file
      --bim                                   Path to bim file
      --fam                                   Path to fam file
   
    Optional:
      --mind  [0.02 (default), 0-1]
      --indep_pairwise ["50 5 0.2" (default), ""]
      --variant_geno [0.02 (default), 0-1]
      --hwe  [1e-7 (default), 0-1]
      --maf [0.05 (default), 0-1]
      --missingness [1e-7(default), 0-1]
	  --sexcheck [true (default),false]
	  --keep_sex_chroms [true (default),false]
	  --pihat [0.125 (default), 0-1]
	  --pca_covars [3 (default), 1-20]

    
  Population stratification options:
    Mandatory:
      --qc

    Optional:
      --maf [0.05 (default), 0-1]
      --indep-pairwise ["50 5 0.2" (default), ""]
      --racefile [super (default), sub]
	  --parfile [false (default), parfile.txt]
	  --racecode [""(default), "EUR"/"AFR"/"SAS"... ]

  Pre-imputation workflow options:
    Mandatory:
      --qc

  Imputation workflow options:
    Mandatory:
      --qc

    Optional:
      --impute_maf [0.01 (default), 0-1]
	  --info [0.7 (default), 0-1]
	
  Post-imputation workflow options:
    Mandatory:
	  --vcf
	  --fam
	
	Optional:
	  --impute_maf [0.01 (default), 0-1]
	  --info [0.7 (default), 0-1]

  GWAS workflow options:
    Mandatory:
      --qc

    Optional:
      --pop_strat
      --impute  
	  --covar_file [false (default), covar.txt]
	  --pca_covars [3 (default), 1-20]
  """.stripIndent()
}
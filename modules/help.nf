def printHelp() {
  log.info"""
  Usage:
    nextflow run main.nf -profile (docker,conda) ( --convertBuild | --qc | --popStrat | --impute | --gwas ) [workflow-options]
  Description:
    Make your single-nucleotide polymorphisms cute with snpQT. A pipeline for quality control, population stratification, imputation, and GWAS of human genomic variants.
  
  Full documentation available at: https://snpqt.readthedocs.org

  Nextflow arguments:
    -profile 
  
  Main workflow arguments:
    --convertBuild                            Run the build conversion workflow
    --qc                                      Run the quality control workflow
    --popStrat                                Run the population stratification workflow
    --impute                                  Run the imputation & post-imputation quality control workflows
    --gwas                                    Run the genome wide association study workflow
    --download_db                             Run the reference database setup workflow

    --results                                 The output directory for results files (default: snpQT/results)

  Quality control workflow options:
    Mandatory:
      --fam
      --bed
      --bim (if input not --vcf)
      --vcf (if input not --bed & --bim)
    Optional:
      --mind  (Default: 0.02)
      --indep_pairwise (Default: "50 5 0.2")
      --variant_geno (Default: 0.02)
      --hwe  (Default: 1e-7)
      --maf (Default: 0.05)
      --missingness (Default: 10e-7)

    Database download workflow options:
      --download-db core                      Download the core reference files to enable --convertBuild, --qc, --popStrat, and --gwas
      --download-db impute                    Download additional reference files to enable --impute

  Population stratification options:
    Mandatory:
      --qc
  
    Optional:
      --maf (default: 0.05)
      --indep-pairwise (Default: "50 5 0.2")
      --racefile ([super (default), sub])

  Imputation and post-imputation options:
    Mandatory:
      --qc

    Optional:
      --info (default: 0.7)

  GWAS options:
    Mandatory:
      --qc
      -popStrat

    Optional:
      --impute  
  """.stripIndent()
}
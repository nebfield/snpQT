def printHelp() {
  log.info"""
  Usage:
    nextflow run main.nf -profile (singularity,docker,conda) ( --buildConversion | --qc | --popStrat | --impute | --postImpute | --gwas ) [workflow-options]
  Description:
    Make your single-nucleotide polymorphisms cute with snpQT. A pipeline for quality control, population stratification, imputation, and GWAS of human genomic variants.
  
  Full documentation available at: https://snpqt.readthedocs.org
  """.stripIndent()
}
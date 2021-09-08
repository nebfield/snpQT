def printHelp() {
  log.info"""
  Usage:
    nextflow run main.nf -profile ([standard,cluster][singularity,docker,conda]) -params-file parameters.yaml
  Description:
    Make your single-nucleotide polymorphisms cute with snpQT. A pipeline for quality control, population stratification, imputation, and GWAS of human genomic variants.
  
  Full documentation with examples, including parameter files, available at: https://snpqt.readthedocs.org

  Nextflow arguments:
        -profile      A comma separated list e.g. standard,singularity or cluster,conda
	-params-file  Path to a YAML file containing all pipeline parameters e.g. parameters.yaml
        -resume       Attempt to continue pipeline from previous checkpoint
  """.stripIndent()
}

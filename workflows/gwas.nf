// GWAS workflow
nextflow.preview.dsl = 2

// import modules
include {run_gwas} from '../modules/gwas.nf' // E1
include {plot} from '../modules/gwas.nf' // E2

workflow gwas {
  take:
    ch_bed
    ch_bim
    ch_fam
    covar
    
  main:
    run_gwas(ch_bed, ch_bim, ch_fam, covar)
    run_gwas.out.logistic.flatten()
      .map { file -> tuple(file.simpleName, file) }
      .set{logistic}
    plot(logistic)
}
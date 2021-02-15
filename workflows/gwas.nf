// GWAS workflow
nextflow.preview.dsl = 2

// import modules
include {run_gwas} from '../modules/gwas.nf' // E1
include {plot} from '../modules/gwas.nf' // E2
include {parse_logs} from '../modules/qc.nf'
include {report} from '../modules/qc.nf'

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
      .set{ logistic }

    plot(logistic)
    parse_logs("gwas", run_gwas.out.log, "gwas_log.txt")
    plot.out.qqplot
      .concat(plot.out.manhattan, parse_logs.out.figure)
      .collect()
      .set{ figures }
    Channel
      .fromPath("$baseDir/bootstrap/gwas_report.Rmd", checkIfExists: true)
      .set{ rmd }
    report("gwas", figures, rmd)  
}

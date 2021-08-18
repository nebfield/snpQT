// GWAS workflow
nextflow.enable.dsl = 2

// import modules
include {run_gwas} from '../modules/gwas.nf' // I1
include {plot} from '../modules/gwas.nf' // I1
include {parse_logs} from '../modules/qc.nf' // E13
include {report} from '../modules/qc.nf' // E14

workflow gwas {
  take:
    ch_bed
    ch_bim
    ch_fam
    covar
    
  main:
    run_gwas(ch_bed, ch_bim, ch_fam, covar)
    run_gwas.out.gwas.flatten()
      .map { file -> tuple(file.simpleName, file) }
      .set{ gwas }
	run_gwas.out.log.flatten()
      .map { file -> tuple(file.simpleName, file) }
      .set{ logs }
	
	gwas
	  .join(logs)
	  .set{ gwas_files }
	  
    plot(gwas_files)
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

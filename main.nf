#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {printHelp} from './modules/help.nf'

// import subworkflows
include {buildConversion} from './workflows/buildConversion.nf'
include {qc} from './workflows/qc.nf'
include {popStrat} from './workflows/popStrat.nf'
include {imputation} from './workflows/imputation.nf'
include {postImputation} from './workflows/postImputation.nf'
include {gwas} from './workflows/gwas.nf'

// initialise default parameters
params.bed = false
params.bim = false
params.fam = false
params.vcf = false

params.convertBuild = false
params.qc = false
params.popStrat = false
params.impute = false
params.postImpute = false
params.gwas = false
params.help = false

// todo: error checking input configuration
if (params.help) {
  printHelp()
  System.exit(0)
}

if (params.convertBuild ) {
  if (!params.vcf) {
    println("Please supply a vcf.gz file for build conversion with --vcf")
  }
  if (params.bed && !params.bim || !params.bed && params.bim ) {
    println("--bed and --bim must be supplied together")
    println("Use --help to print help")
    System.exit(1)
  }
} else if (params.qc) {
  if (params.qc && !params.convertBuild) {
    if(params.vcf) {
      println("QC module is not compatible with direct VCF input") // TODO: fix this 
      println("Please supply bim, bed, and fam files with --bim, --bed, and --fam")
      println("Use --help to print help")
      System.exit(1)
    }
  }
} else if (params.gwas) {
    if(!params.qc || !params.popStrat) {
      println("GWAS module requires qc and popStrat")
      println("Please rerun with --qc and --popStrat")
      println("Use --help to print help")
      System.exit(1)
    }
}

// main workflow
workflow {
  Channel
    .fromPath(params.db, checkIfExists: true)
    .set{ ch_db }

  if ( params.convertBuild ) {
    Channel
      .fromPath(params.vcf, checkIfExists: true)
      .set{ ch_vcf }
  }
  
  if ( !params.convertBuild && params.qc ) {
    Channel
      .fromPath(params.bed, checkIfExists: true)
      .set{ ch_bed }
    Channel
      .fromPath(params.bim, checkIfExists: true)
      .set{ ch_bim }
    Channel
      .fromPath(params.fam, checkIfExists: true)
      .set{ ch_fam }
  } else if (params.convertBuild && params.qc) {
    Channel
      .fromPath(params.fam, checkIfExists: true)
      .set{ ch_fam }
  }

  if ( params.popStrat ) {
    Channel
      .fromPath(params.db + '/all_phase3_10.bed', checkIfExists: true )
      .set{ ch_ref_bed }
    Channel
      .fromPath(params.db + '/all_phase3_10.bim', checkIfExists: true )
      .set{ ch_ref_bim }
    Channel
      .fromPath(params.db + '/all_phase3_10.fam', checkIfExists: true)
      .set{ ch_ref_fam } 
  }

  if ( params.impute ) {
  }

  main:
    if ( params.convertBuild && !params.qc) {
      buildConversion(ch_vcf, ch_db)
    }
    
    if ( params.convertBuild && params.qc ) {
      buildConversion(ch_vcf, ch_db)
      qc(buildConversion.out.bed, buildConversion.out.bim, ch_fam, ch_db)
    }

    if (params.qc && params.popStrat) {
      popStrat(qc.out.bed, qc.out.bim, qc.out.fam, ch_ref_bed, ch_ref_bim, ch_ref_fam, ch_db)
    }

    if (params.qc && params.impute) {
      imputation(qc.out.bed, qc.out.bim, qc.out.fam, ch_db)
    }

    if (params.impute && params.postImpute && params.qc) {
      postImputation(imputation.out.imputed, qc.out.fam)
    }

    if (params.qc && params.popStrat && params.gwas) {
      gwas(qc.out.bed, qc.out.bim, qc.out.fam, popStrat.out.covar)
    }
    
}

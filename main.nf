#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {printHelp} from './modules/help.nf'

// import subworkflows
include {buildConversion} from './workflows/buildConversion.nf'
include {sample_qc} from './workflows/sample_qc.nf'
include {variant_qc} from './workflows/variant_qc.nf'
include {popStrat} from './workflows/popStrat.nf'
include {imputation} from './workflows/imputation.nf'
include {postImputation} from './workflows/postImputation.nf'
include {gwas} from './workflows/gwas.nf'

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
  // set up input channels

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
  
  main:
    if ( params.convertBuild) {
      buildConversion(ch_vcf)
    }
    
    if ( params.convertBuild && params.qc && !params.popStrat) {
      sample_qc(buildConversion.out.bed, buildConversion.out.bim, ch_fam)
      variant_qc(sample_qc.out.bed, sample_qc.out.bim, sample_qc.out.fam)
    } else if ( !params.convertBuild && params.qc ) {
      sample_qc(ch_ref_bed, ch_ref_bim, ch_ref_fam)
    }

    if (params.qc && params.popStrat) {
      sample_qc(buildConversion.out.bed, buildConversion.out.bim, ch_fam)
      popStrat(sample_qc.out.bed, sample_qc.out.bim, sample_qc.out.fam)
      variant_qc(popStrat.out.bed, popStrat.out.bim, popStrat.out.fam)      
    }

    if (params.qc && params.impute) {
      // imputation(qc.out.bed, qc.out.bim, qc.out.fam, ch_db)
    }

    if (params.impute && params.postImpute && params.qc) {
      // postImputation(imputation.out.imputed, qc.out.fam)
    }

    if (params.qc && params.popStrat && params.gwas) {
      // gwas(qc.out.bed, qc.out.bim, qc.out.fam, popStrat.out.covar)
    }
    
}

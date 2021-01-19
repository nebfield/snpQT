#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import subworkflows
include {buildConversion} from './workflows/buildConversion.nf'
include {qc} from './workflows/qc.nf'

// initialise default parameters
params.bed = false
params.bim = false
params.fam = false
params.vcf = false

params.convertBuild = false
params.qc = false
params.help = true

// todo: error checking input configuration
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
      .set{ ch_inbed }
    Channel
      .fromPath(params.bim, checkIfExists: true)
      .set{ ch_inbim }
    Channel
      .fromPath(params.fam, checkIfExists: true)
      .set{ ch_infam }
  } else if (params.convertBuild && params.qc) {
    Channel
      .fromPath(params.fam, checkIfExists: true)
      .set{ ch_infam }
  }

  main:
    if ( params.convertBuild && !params.qc) {
      buildConversion(ch_vcf, ch_db)
    }
    
    if ( params.convertBuild && params.qc ) {
      buildConversion(ch_vcf, ch_db)
      qc(buildConversion.out.bed, buildConversion.out.bim, ch_infam, ch_db)
    }
}

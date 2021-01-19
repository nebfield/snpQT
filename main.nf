#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import subworkflows
include {qc} from './workflows/qc.nf'

// todo: error checking parameters

// main workflow
workflow {
  Channel
    .fromPath(params.db)
    .set{ ch_db }
  if ( params.qc ) {
    Channel
      .fromPath(params.inbed)
      .set{ ch_inbed }
    Channel
      .fromPath(params.inbim)
      .set{ ch_inbim }
    Channel
      .fromPath(params.infam)
      .set{ ch_infam }
  }

  main:
    if ( params.qc ) {
      qc(ch_inbed, ch_inbim, ch_infam, ch_db)
    }
}
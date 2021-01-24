// Variant quality control workflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {mpv} from '../modules/qc.nf' // B8
include {plot_mpv} from '../modules/qc.nf' // B8
include {hardy} from '../modules/qc.nf' // B9
include {plot_hardy} from '../modules/qc.nf' // B9
include {maf} from '../modules/qc.nf' // B10
include {plot_maf} from '../modules/qc.nf' // B10
include {test_missing} from '../modules/qc.nf' // B11
include {plot_missing_by_cohort} from '../modules/qc.nf' // B11

// workflow component for snpqt pipeline
workflow variant_qc {
  take:
    ch_inbed
    ch_inbim
    ch_infam

  main:
    mpv(ch_inbed, ch_inbim, ch_infam)
    plot_mpv(mpv.out.lmiss)
    hardy(mpv.out.bed, mpv.out.bim, mpv.out.fam)
    plot_hardy(hardy.out.sub, hardy.out.zoom)
    maf(hardy.out.bed, hardy.out.bim, hardy.out.fam)
    plot_maf(maf.out.frq)
    test_missing(maf.out.bed, maf.out.bim, maf.out.fam)
    plot_missing_by_cohort(test_missing.out.missing)

  emit:
    bed = test_missing.out.bed
    bim = test_missing.out.bim
    fam = test_missing.out.fam
} 
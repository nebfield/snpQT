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
include {parse_logs} from '../modules/qc.nf'
include {pca_covariates} from '../modules/popStrat.nf' // C10

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
    Channel
      .fromPath("$baseDir/db/PCA.exclude.regions.b37.txt", checkIfExists: true)
      .set{ exclude }
    pca_covariates(test_missing.out.bed, test_missing.out.bim, test_missing.out.fam, exclude)
    logs = mpv.out.log.concat(hardy.out.log, maf.out.log, test_missing.out.log).collect()
    parse_logs("qc", logs, "variant_qc.log")
 
  emit:
    bed = test_missing.out.bed
    bim = test_missing.out.bim
    fam = test_missing.out.fam
    covar = pca_covariates.out.covar
} 

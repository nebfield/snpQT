// Variant quality control workflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {mpv} from '../modules/qc.nf' // B8
include {plot_mpv as mpv_before; plot_mpv as mpv_after} from '../modules/qc.nf' // B8
include {hardy} from '../modules/qc.nf' // B9
include {plot_hardy as ph_before; plot_hardy as ph_after} from '../modules/qc.nf' // B9
include {maf} from '../modules/qc.nf' // B10
include {plot_maf as pm_before; plot_maf as pm_after} from '../modules/qc.nf' // B10
include {test_missing} from '../modules/qc.nf' // B11
include {plot_missing_by_cohort as pmc_before; plot_missing_by_cohort as pmc_after} from '../modules/qc.nf' // B11
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
    mpv_before(mpv.out.lmiss_before, params.variant_geno, "before")
    mpv_after(mpv.out.lmiss_after, params.variant_geno, "after")
    hardy(mpv.out.bed, mpv.out.bim, mpv.out.fam)
    ph_before(hardy.out.sub_before, hardy.out.zoom_before, params.hwe, "before")
    ph_after(hardy.out.sub_after, hardy.out.zoom_after, params.hwe, "after")
    maf(hardy.out.bed, hardy.out.bim, hardy.out.fam)
    pm_before(maf.out.before, params.maf, "before")
    pm_after(maf.out.after, params.maf, "after")
    test_missing(maf.out.bed, maf.out.bim, maf.out.fam)
    pmc_before(test_missing.out.before, params.missingness, "before")
    pmc_after(test_missing.out.after, params.missingness, "after")
   
   Channel
      .fromPath("$baseDir/db/PCA.exclude.regions.b37.txt", checkIfExists: true)
      .set{ exclude }
    pca_covariates(test_missing.out.bed, test_missing.out.bim, test_missing.out.fam, exclude)
    logs = mpv.out.log.concat(hardy.out.log, maf.out.log, test_missing.out.log).collect()
    parse_logs("qc", logs, "variant_qc_log.txt")
 
  emit:
    bed = test_missing.out.bed
    bim = test_missing.out.bim
    fam = test_missing.out.fam
    covar = pca_covariates.out.covar
} 

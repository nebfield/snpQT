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
include {pca} from '../modules/qc.nf' // B12
include {pca_covariates} from '../modules/qc.nf' // B13
include {plot_pca_user_data} from '../modules/qc.nf' // B14
include {report} from '../modules/qc.nf'

// workflow component for snpqt pipeline
workflow variant_qc {
  take:
    ch_inbed
    ch_inbim
    ch_infam

  main:
    mpv(ch_inbed, ch_inbim, ch_infam)
    plot_mpv(mpv.out.lmiss_before, mpv.out.lmiss_after, params.variant_geno)
    hardy(mpv.out.bed, mpv.out.bim, mpv.out.fam)
    plot_hardy(hardy.out.sub_before, hardy.out.zoom_before, hardy.out.sub_after, hardy.out.zoom_after, params.hwe)
    maf(hardy.out.bed, hardy.out.bim, hardy.out.fam)
    plot_maf(maf.out.before, maf.out.after, params.maf)
    test_missing(maf.out.bed, maf.out.bim, maf.out.fam)
    plot_missing_by_cohort(test_missing.out.before, test_missing.out.after, params.missingness)
    Channel
      .fromPath("$baseDir/db/PCA.exclude.regions.b37.txt", checkIfExists: true)
      .set{ exclude }
    pca(test_missing.out.bed, test_missing.out.bim, test_missing.out.fam, exclude)
    plot_pca_user_data(pca.out.eigenvec_user, pca.out.fam)
	if (params.covar_file == false) {
	  pca_covariates(pca.out.eigenvec_user)
	  covar = pca_covariates.out.covar
    }else{
	  Channel
        .fromPath(params.covar_file, checkIfExists: true)
        .set{ covar }
    }
	logs = mpv.out.log.concat(hardy.out.log, maf.out.log, test_missing.out.log).collect()
    parse_logs("qc", logs, "variant_qc_log.txt")
    figures = plot_mpv.out.figure
      .concat(plot_hardy.out.figure, plot_maf.out.figure, plot_missing_by_cohort.out.figure, plot_pca_user_data.out.figure, plot_pca_user_data.out.rds ,parse_logs.out.figure)
      .collect()
    Channel
      .fromPath("$baseDir/bootstrap/variant_report.Rmd", checkIfExists: true)
      .set{ rmd }
    report("qc", figures, rmd)
 
  emit:
    bed = test_missing.out.bed
    bim = test_missing.out.bim
    fam = test_missing.out.fam
    covar
} 

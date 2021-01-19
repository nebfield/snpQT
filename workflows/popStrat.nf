// Population stratification workflow
nextflow.preview.dsl = 2

// import modules
include {filter_maf} from '../modules/popStrat.nf' // C3
include {run_snpflip} from '../modules/popStrat.nf' // C4 
include {flip_snps} from '../modules/popStrat.nf' // C4
include {align} from '../modules/popStrat.nf' // C5
include {merge} from '../modules/popStrat.nf' // C6
include {pca_prep} from '../modules/popStrat.nf' // C7
include {racefile} from '../modules/popStrat.nf' // C7
include {eigensoft} from '../modules/popStrat.nf' // C8
include {plot_pca} from '../modules/popStrat.nf' // C8
include {extract_homogenous} from '../modules/popStrat.nf' // C9
include {pca_covariates} from '../modules/popStrat.nf' // C10

// workflow component for snpqt pipeline
workflow popStrat {
  take:
    ch_bed
    ch_bim
    ch_fam
    ref_bed
    ref_bim
    ref_fam
    ch_db

  main:
    filter_maf(ch_bed, ch_bim, ch_fam, ch_db)
    run_snpflip(filter_maf.out.bed, filter_maf.out.bim, filter_maf.out.fam, ch_db)
    flip_snps(filter_maf.out.bed, filter_maf.out.bim, filter_maf.out.fam, run_snpflip.out.rev, run_snpflip.out.ambig)
    align(flip_snps.out.bed, flip_snps.out.bim, flip_snps.out.fam, ref_bed, ref_bim, ref_fam)
    merge(align.out.bed, align.out.bim, align.out.fam, ref_bed, ref_bim, ref_fam)
    pca_prep(merge.out.bed, merge.out.bim, merge.out.fam, ch_db)
    racefile(ch_db)
    eigensoft(pca_prep.out.bed, pca_prep.out.bim, pca_prep.out.fam, racefile.out.super, filter_maf.out.fam)
    plot_pca(eigensoft.out.eigenvec, eigensoft.out.merged_racefile)
    extract_homogenous(filter_maf.out.bed, filter_maf.out.bim, filter_maf.out.fam, eigensoft.out.keep_samples)
    pca_covariates(extract_homogenous.out.bed, extract_homogenous.out.bim, extract_homogenous.out.fam, ch_db)

  emit:
    covar = pca_covariates.out.covar
} 
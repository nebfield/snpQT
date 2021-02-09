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
include {pca_plink} from '../modules/popStrat.nf' // C8 
include {extract_homogenous} from '../modules/popStrat.nf' // C9
include {parse_logs} from '../modules/qc.nf'

// workflow component for snpqt pipeline
workflow popStrat {
  take:
    ch_bed
    ch_bim
    ch_fam

  main:
    Channel
      .fromPath("$baseDir/db/all_phase3_10.bed", checkIfExists: true)
      .set{ ref_bed }
    Channel
      .fromPath("$baseDir/db/all_phase3_10.bim", checkIfExists: true)
      .set{ ref_bim }
    Channel
      .fromPath("$baseDir/db/all_phase3_10.fam", checkIfExists: true)
      .set{ ref_fam } 
    Channel
      .fromPath("$baseDir/db/PCA.exclude.regions.b37.txt", checkIfExists: true)
      .set{ exclude }
    filter_maf(ch_bed, ch_bim, ch_fam, exclude)
    Channel
      .fromPath("$baseDir/db/h37_squeezed.fasta", checkIfExists: true)
      .set { g37 }
    run_snpflip(filter_maf.out.bed, filter_maf.out.bim, filter_maf.out.fam, g37)
    flip_snps(filter_maf.out.bed, filter_maf.out.bim, filter_maf.out.fam, run_snpflip.out.rev, run_snpflip.out.ambig)
    align(flip_snps.out.bed, flip_snps.out.bim, flip_snps.out.fam, ref_bed, ref_bim, ref_fam)
    merge(align.out.bed, align.out.bim, align.out.fam, ref_bed, ref_bim, ref_fam)
    pca_prep(merge.out.bed, merge.out.bim, merge.out.fam, exclude)
    Channel
      .fromPath("$baseDir/db/integrated_call_samples_v3.20130502.ALL.panel", checkIfExists: true)
      .set{ panel }
    racefile(panel)
    if (params.racefile == "super") {
      rf = racefile.out.super
    } else if (params.racefile == "sub") {
      rf = racefile.out.sub
    }
    eigensoft(pca_prep.out.bed, pca_prep.out.bim, pca_prep.out.fam, rf, filter_maf.out.fam)
    plot_pca(eigensoft.out.eigenvec, eigensoft.out.merged_racefile)
    pca_plink(pca_prep.out.bed, pca_prep.out.bim, pca_prep.out.fam, eigensoft.out.eigenvec)
    extract_homogenous(ch_bed, ch_bim, ch_fam, eigensoft.out.keep_samples)
    logs = filter_maf.out.log.concat(flip_snps.out.log, align.out.log, merge.out.log, pca_prep.out.log, extract_homogenous.out.log).collect()
    parse_logs("popStrat", logs, "popStrat_log.txt")

  emit:
    bed = extract_homogenous.out.bed
    bim = extract_homogenous.out.bim
    fam = extract_homogenous.out.fam
} 

// Quality control workflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {variant_missingness} from '../modules/qc.nf' // B1
include {individual_missingness} from '../modules/qc.nf' // B2
include {plot_missingness} from '../modules/qc.nf' // B2
include {check_sex} from '../modules/qc.nf' // B3
include {plot_sex} from '../modules/qc.nf' // B3
include {extract_autosomal} from '../modules/qc.nf' // B4
include {heterozygosity_rate} from '../modules/qc.nf' // B5
include {plot_heterozygosity} from '../modules/qc.nf' // B5
include {heterozygosity_prune} from '../modules/qc.nf' // B5
include {relatedness} from '../modules/qc.nf' // B6
include {missing_phenotype} from '../modules/qc.nf' // B7
include {mpv} from '../modules/qc.nf' // B8
include {plot_mpv} from '../modules/qc.nf' // B8
include {hardy} from '../modules/qc.nf' // B9
include {plot_hardy} from '../modules/qc.nf' // B9
include {maf} from '../modules/qc.nf' // B10
include {plot_maf} from '../modules/qc.nf' // B10
include {test_missing} from '../modules/qc.nf' // B11

// workflow component for snpqt pipeline
workflow qc {
  take:
    ch_inbed
    ch_inbim
    ch_infam
    ch_db

  main:
    variant_missingness(ch_inbed, ch_inbim, ch_infam)
    individual_missingness(variant_missingness.out.bed, variant_missingness.out.bim, variant_missingness.out.fam)
    plot_missingness(individual_missingness.out.imiss)
    check_sex(individual_missingness.out.bed, individual_missingness.out.bim, individual_missingness.out.fam)
    plot_sex(check_sex.out.sexcheck)
    extract_autosomal(check_sex.out.bed, check_sex.out.bim, check_sex.out.fam)
    heterozygosity_rate(extract_autosomal.out.bed, extract_autosomal.out.bim, extract_autosomal.out.fam, ch_db)
    plot_heterozygosity(heterozygosity_rate.out.het)
    heterozygosity_prune(extract_autosomal.out.bed, extract_autosomal.out.bim, extract_autosomal.out.fam, plot_heterozygosity.out.failed)
    relatedness(heterozygosity_prune.out.bed, heterozygosity_prune.out.bim, heterozygosity_prune.out.fam, heterozygosity_rate.out.ind_snps, individual_missingness.out.imiss)
    missing_phenotype(relatedness.out.bed, relatedness.out.bim, relatedness.out.fam)
    mpv(missing_phenotype.out.bed, missing_phenotype.out.bim, missing_phenotype.out.fam)
    plot_mpv(mpv.out.lmiss)
    hardy(mpv.out.bed, mpv.out.bim, mpv.out.fam)
    plot_hardy(hardy.out.sub, hardy.out.zoom)
    maf(hardy.out.bed, hardy.out.bim, hardy.out.fam)
    plot_maf(maf.out.frq)
    test_missing(maf.out.bed, maf.out.bim, maf.out.fam)

  emit:
    qc_bed = test_missing.out.bed
    qc_bim = test_missing.out.bim
    qc_fam = test_missing.out.fam
    reports = combine(plot_missingness.out.figure, plot_sex.out.figure, plot_heterozygosity.out.figure, plot_mpv.out.figure, plot_hardy.out.figure, plot_maf.out.figure)
} 
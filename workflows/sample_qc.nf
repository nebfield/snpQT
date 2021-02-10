// Quality control workflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {variant_missingness} from '../modules/qc.nf' // B1
include {individual_missingness} from '../modules/qc.nf' // B2
include {plot_missingness as pm_before; plot_missingness as pm_after} from '../modules/qc.nf' // B2
include {check_sex} from '../modules/qc.nf' // B3
include {plot_sex} from '../modules/qc.nf' // B3
include {extract_autosomal} from '../modules/qc.nf' // B4
include {heterozygosity_rate} from '../modules/qc.nf' // B5
include {plot_heterozygosity as ph_before; plot_heterozygosity as ph_after} from '../modules/qc.nf' // B5
include {heterozygosity_prune} from '../modules/qc.nf' // B5
include {relatedness} from '../modules/qc.nf' // B6
include {missing_phenotype} from '../modules/qc.nf' // B7
include {parse_logs} from '../modules/qc.nf'

// workflow component for snpqt pipeline
workflow sample_qc {
  take:
    ch_inbed
    ch_inbim
    ch_infam
    
  main:
    variant_missingness(ch_inbed, ch_inbim, ch_infam)
    individual_missingness(variant_missingness.out.bed, variant_missingness.out.bim, variant_missingness.out.fam)
    // plot missingness
    pm_before(individual_missingness.out.imiss_before, params.mind, "before")
    pm_after(individual_missingness.out.imiss_after, params.mind, "after")    
    Channel
      .fromPath("$baseDir/db/PCA.exclude.regions.b37.txt", checkIfExists: true)
      .set{exclude}      
    if (params.sexcheck) {
      check_sex(individual_missingness.out.bed, individual_missingness.out.bim, individual_missingness.out.fam)
      plot_sex(check_sex.out.sexcheck)
      extract_autosomal(check_sex.out.bed, check_sex.out.bim, check_sex.out.fam)
      heterozygosity_rate(extract_autosomal.out.bed, extract_autosomal.out.bim, extract_autosomal.out.fam, exclude)
      // plot heterozygosity 
      ph_before(heterozygosity_rate.out.het,"before")   
      heterozygosity_prune(extract_autosomal.out.bed, extract_autosomal.out.bim, extract_autosomal.out.fam, ph_before.out.failed, heterozygosity_rate.out.ind_snps)
      ph_after(heterozygosity_prune.out.het,"after")
	  } else {
      heterozygosity_rate(individual_missingness.out.bed, individual_missingness.out.bim, individual_missingness.out.fam, exclude)
      // plot heterozygosity 
      ph_before(heterozygosity_rate.out.het,"before")   
      heterozygosity_prune(individual_missingness.out.bed, individual_missingness.out.bim, individual_missingness.out.fam, ph_before.out.failed)
      ph_after(heterozygosity_prune.out.het,"after")
    }
    relatedness(heterozygosity_prune.out.bed, heterozygosity_prune.out.bim, heterozygosity_prune.out.fam, heterozygosity_rate.out.ind_snps, individual_missingness.out.imiss_after)
    missing_phenotype(relatedness.out.bed, relatedness.out.bim, relatedness.out.fam)
    logs = variant_missingness.out.log.concat(individual_missingness.out.log, check_sex.out.log, extract_autosomal.out.log, heterozygosity_prune.out.log, relatedness.out.log, missing_phenotype.out.log).collect()
    parse_logs("qc", logs, "sample_qc_log.txt")
  

  emit:
    bed = missing_phenotype.out.bed
    bim = missing_phenotype.out.bim
    fam = missing_phenotype.out.fam
} 

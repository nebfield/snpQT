// Post-imputation QC workflow
nextflow.preview.dsl = 2

// import modules
include {merge_imp} from '../modules/postImputation.nf' // E1
include {filter_imp} from '../modules/postImputation.nf' // E2
include {duplicates_cat1} from '../modules/postImputation.nf' // E3
include {duplicates_cat2} from '../modules/postImputation.nf' // E3
include {duplicates_cat3} from '../modules/postImputation.nf' // E3
include {update_phenotype} from '../modules/postImputation.nf' // E4
include {parse_logs} from '../modules/qc.nf'

workflow postImputation {
  take:
    ch_imp
    ch_fam
    
  main:
    merge_imp(ch_imp)
    filter_imp(merge_imp.out.vcf)
    duplicates_cat1(filter_imp.out.bed, filter_imp.out.bim, filter_imp.out.fam)
    duplicates_cat2(duplicates_cat1.out.bed, duplicates_cat1.out.bim, duplicates_cat1.out.fam)
    duplicates_cat3(duplicates_cat2.out.bed, duplicates_cat2.out.bim, duplicates_cat2.out.fam)
    update_phenotype(duplicates_cat3.out.bed, duplicates_cat3.out.bim, duplicates_cat3.out.fam, ch_fam)
    //logs = filter_imp.out.log.concat(duplicates_cat1.out.log, duplicates_cat2.log, duplicates_cat3.log, update_phenotype.out.log).collect()
    //parse_logs("imputation", logs, "postImpute_log.txt")

  emit:
    bed = update_phenotype.out.bed
    bim = update_phenotype.out.bim
    fam = update_phenotype.out.fam
}

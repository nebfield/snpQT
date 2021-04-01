// Post-imputation QC workflow
nextflow.preview.dsl = 2

// import modules
include {filter_imp} from '../modules/postImputation.nf' // H1
include {filter_maf} from '../modules/postImputation.nf' // H2
include {duplicates_cat1} from '../modules/postImputation.nf' // H3
include {duplicates_cat2} from '../modules/postImputation.nf' // H4
include {duplicates_cat3} from '../modules/postImputation.nf' // H5
include {update_ids} from '../modules/postImputation.nf' // H6
include {update_phenotype} from '../modules/postImputation.nf' // H7
include {parse_logs} from '../modules/qc.nf' // E13

workflow postImputation {
  take:
    ch_imp
    ch_fam
    
  main:
    filter_imp(ch_imp)
    filter_maf(filter_imp.out.bed,filter_imp.out.bim,filter_imp.out.fam)
    duplicates_cat1(filter_maf.out.bed, filter_maf.out.bim, filter_maf.out.fam)
    duplicates_cat2(duplicates_cat1.out.bed, duplicates_cat1.out.bim, duplicates_cat1.out.fam)
    duplicates_cat3(duplicates_cat2.out.bed, duplicates_cat2.out.bim, duplicates_cat2.out.fam)
    update_ids(duplicates_cat3.out.bed, duplicates_cat3.out.bim, duplicates_cat3.out.fam, ch_fam)
	update_phenotype(update_ids.out.bed, update_ids.out.bim, update_ids.out.fam, ch_fam)
    logs = filter_imp.out.log.concat(filter_maf.out.log, duplicates_cat1.out.log, duplicates_cat2.out.log, duplicates_cat3.out.log, update_ids.out.log, update_phenotype.out.log).collect()
    parse_logs("post_imputation", logs, "post_impute_log.txt")

  emit:
    bed = update_phenotype.out.bed
    bim = update_phenotype.out.bim
    fam = update_phenotype.out.fam
}

// Post-imputation QC workflow
nextflow.preview.dsl = 2

// import modules
include {merge_imp} from '../modules/postImputation.nf' // E1
include {filter_imp} from '../modules/postImputation.nf' // E2
include {filter_maf} from '../modules/postImputation.nf' // E3
include {duplicates_cat1} from '../modules/postImputation.nf' // E4
include {duplicates_cat2} from '../modules/postImputation.nf' // E5
include {duplicates_cat3} from '../modules/postImputation.nf' // E6
include {update_phenotype} from '../modules/postImputation.nf' // E7
include {parse_logs} from '../modules/qc.nf'

workflow postImputation {
  take:
    ch_imp
    ch_fam
    
  main:
    merge_imp(ch_imp)
	annotate_missing(merge_imp.out.vcf)
    filter_imp(annotate_missing.out.bed,annotate_missing.out.bim,annotate_missing.out.fam)
	filter_maf(filter_imp.out.bed,filter_imp.out.bim,filter_imp.out.fam)
    duplicates_cat1(filter_maf.out.bed, filter_maf.out.bim, filter_maf.out.fam)
    duplicates_cat2(duplicates_cat1.out.bed, duplicates_cat1.out.bim, duplicates_cat1.out.fam)
    duplicates_cat3(duplicates_cat2.out.bed, duplicates_cat2.out.bim, duplicates_cat2.out.fam)
    update_phenotype(duplicates_cat3.out.bed, duplicates_cat3.out.bim, duplicates_cat3.out.fam, ch_fam)
    logs = annotate_missing.out.log.concat(filter_imp.out.log,filter_maf.out.log, duplicates_cat1.out.log, duplicates_cat2.out.log, duplicates_cat3.out.log, update_phenotype.out.log).collect()
    parse_logs("imputation", logs, "postImpute_log.txt")

  emit:
    bed = update_phenotype.out.bed
    bim = update_phenotype.out.bim
    fam = update_phenotype.out.fam
}

// Build conversion workflow
nextflow.preview.dsl = 2

// import modules
include {dictionary} from '../modules/buildConversion.nf' // A4
include {num_to_chr} from '../modules/buildConversion.nf' // A6
include {liftover} from '../modules/buildConversion.nf' // A7
include {chr_to_num} from '../modules/buildConversion.nf' // A8 & A9
include {vcf_to_plink} from '../modules/buildConversion.nf' // A10

// workflow component for snpqt pipeline
workflow buildConversion {
  take:
    ch_vcf
    ch_db

  main:
    dictionary(ch_db)
    num_to_chr(ch_vcf, ch_db)
    liftover(num_to_chr.out.vcf, ch_db, dictionary.out.dict)
    chr_to_num(liftover.out.vcf, ch_db)
    vcf_to_plink(chr_to_num.out.vcf)

  emit:
    bed = vcf_to_plink.out.bed
    bim = vcf_to_plink.out.bim
    fam = vcf_to_plink.out.fam
    
} 
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

  main:
    Channel
      .fromPath("$baseDir/db/picard.jar", checkIfExists: true)
      .set{picard}
    Channel
      .fromPath("$baseDir/db/hg19.fa.gz", checkIfExists: true)
      .set{hg19}
    dictionary(picard, hg19)
    Channel
      .fromPath("$baseDir/db/1toChr1.txt", checkIfExists: true)
      .set{ chr_map }
    num_to_chr(ch_vcf, chr_map)
    Channel
      .fromPath("$baseDir/db/hg38ToHg19.over.chain", checkIfExists: true)
      .set{chain}
    liftover(num_to_chr.out.vcf, picard, hg19, chain, dictionary.out.dict)
    chr_to_num(liftover.out.vcf, chr_map)
    vcf_to_plink(chr_to_num.out.vcf)

  emit:
    bed = vcf_to_plink.out.bed
    bim = vcf_to_plink.out.bim
    fam = vcf_to_plink.out.fam
    
} 
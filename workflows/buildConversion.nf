// Build conversion workflow
nextflow.enable.dsl = 2

// import modules
include {decompress; decompress as decompress_otherBuild} from '../modules/download_db.nf' // A1
include {dictionary} from '../modules/buildConversion.nf' // B1
include {num_to_chr} from '../modules/buildConversion.nf' // B2
include {liftover} from '../modules/buildConversion.nf' // B3
include {chr_to_num} from '../modules/buildConversion.nf' // B4
include {vcf_to_plink} from '../modules/buildConversion.nf' // B5

// workflow component for snpqt pipeline
workflow buildConversion {
  take:
    ch_vcf
	ch_fam

  main:
  
  if (params.input_build == 38 && params.output_build == 37){
    Channel
      .fromPath("${params.db}/hg19.fa.gz", checkIfExists: true)
      .set{hg19}
	decompress(hg19)
    dictionary(decompress.out.file)
    Channel
      .fromPath("${params.db}/1toChr1.txt", checkIfExists: true)
      .set{chr_map}
	num_to_chr(ch_vcf, chr_map)
    Channel
      .fromPath("${params.db}/hg38ToHg19.over.chain", checkIfExists: true)
      .set{chain}
    liftover(num_to_chr.out.vcf, decompress.out.file, chain, dictionary.out.dict)
    chr_to_num(liftover.out.vcf, chr_map)
	}else if (params.input_build == 37 && params.output_build == 38){
		Channel
		  .fromPath("${params.db}/hg38.fa.gz", checkIfExists: true)
		  .set{hg38}
		decompress(hg38)
		dictionary(decompress.out.file)
		Channel
		  .fromPath("${params.db}/1toChr1.txt", checkIfExists: true)
		  .set{chr_map}
		num_to_chr(ch_vcf, chr_map)
		Channel
		  .fromPath("${params.db}/hg19ToHg38.over.chain", checkIfExists: true)
		  .set{chain}
		liftover(num_to_chr.out.vcf, decompress.out.file, chain, dictionary.out.dict)
		chr_to_num(liftover.out.vcf, chr_map)
	}
	
	
    vcf_to_plink(chr_to_num.out.vcf, ch_fam)

  emit:
    bed = vcf_to_plink.out.bed
    bim = vcf_to_plink.out.bim
    fam = vcf_to_plink.out.fam
    
} 
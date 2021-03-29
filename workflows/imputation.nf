// Imputation workflow
nextflow.preview.dsl = 2

// import modules
include {split_user_chrom} from '../modules/imputation.nf' // D12 - D13
include {phasing} from '../modules/imputation.nf' // D14 
include {bcftools_index_chr} from '../modules/imputation.nf' // D15
include {convert_imp5} from '../modules/imputation.nf' // D17 (D16 in container)
include {impute5} from '../modules/imputation.nf' // D18
include {merge_imp} from '../modules/imputation.nf' // D19
include {parse_logs} from '../modules/qc.nf'

// workflow component for snpqt pipeline
workflow imputation {
  take:
    ch_vcf
    ch_idx

  main:
    // OK, things start to get a bit complicated here
    // because we're working with individual chromosomes it's easiest to work with
    // tuples like [chr_num, chr.fa.gz, chr.fa.gz.csi, etc...]
    // that way we can join and combine files by chr_num as needed
    Channel.from(1..22).set{ chrom }
    split_user_chrom(ch_vcf, ch_idx, chrom)
    // maps for shapeit4
    Channel
      .fromPath("$baseDir/db/impute/genetic_maps.b37.tar.gz", checkIfExists: true)
      .set { ch_map }
    user_chroms = split_user_chrom.out.chrom.combine(ch_map)
    phasing(user_chroms)
    phased = bcftools_index_chr(phasing.out.chrom)
    // thousand genome reference data
    Channel
      .fromPath("$baseDir/db/impute/chr*.vcf.gz", checkIfExists: true)
      .map{ f -> [f.baseName.find(/(?<=chr)[A-Z0-9]+/), f] } // regex: chrX, chr11 -> X, 11 
      .set{ thousand_genomes }
    Channel
      .fromPath("$baseDir/db/impute/chr*.vcf.gz.csi", checkIfExists: true)
      .map{ f -> [f.baseName.find(/(?<=chr)[A-Z0-9]+/), f] } // regex: chrX, chr11 -> X, 11 
      .join( thousand_genomes )
      .set{ thousand_genomes_with_idx }
    convert_imp5(thousand_genomes_with_idx)
    impute5(convert_imp5.out.chrom.join(phased).combine(ch_map))
	merge_imp(impute5.out.imputed.collect())
    // logs = fix_duplicates.out.log.concat(to_bcf.out.log).collect()
    // parse_logs("imputation", logs, "imputation_log.txt")

  emit:
    imputed_vcf = merge_imp.out.vcf

} 

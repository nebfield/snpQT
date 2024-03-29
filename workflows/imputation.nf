// Imputation workflow
nextflow.enable.dsl = 2

// import modules
include {annotate_ids} from '../modules/download_db.nf' // A6
include {index_vcf} from '../modules/download_db.nf' // A7
include {split_user_chrom} from '../modules/imputation.nf' // G1
include {phasing} from '../modules/imputation.nf' // G2 
include {bcftools_index_chr} from '../modules/imputation.nf' // G3
include {tabix_chr} from '../modules/imputation.nf' // G4
include {convert_imp5} from '../modules/imputation.nf' // G5
include {impute5} from '../modules/imputation.nf' // G6
include {merge_imp} from '../modules/imputation.nf' // G7
include {parse_logs} from '../modules/qc.nf' // E13

// workflow component for snpqt pipeline
workflow imputation {
  take:
    ch_vcf

  main:
    // annotate id needs a tuple of [id, path], so use a dummy value
    in_vcf = Channel.from('user').combine(ch_vcf)
    annotate_ids(in_vcf)
    index_vcf(annotate_ids.out.vcf)
    // OK, things start to get a bit complicated here
    // because we're working with individual chromosomes it's easiest to work with
    // tuples like [chr_num, chr.fa.gz, chr.fa.gz.tbi, etc...]
    // that way we can join and combine files by chr_num as needed
    Channel.from(1..22).concat(Channel.from('X')).set{ chrom }
    split_user_chrom(annotate_ids.out.vcf, index_vcf.out.idx, chrom)
    // maps for shapeit4
    Channel
      .fromPath("${params.db}/impute/genetic_maps.b37.tar.gz", checkIfExists: true)
      .set { ch_map }
    user_chroms = split_user_chrom.out.chrom.combine(ch_map)
    phasing(user_chroms)
    phased = bcftools_index_chr(phasing.out.chrom)
    // thousand genome reference data
    Channel
      .fromPath("${params.db}/impute/chr*.vcf.gz", checkIfExists: true)
      .map{ f -> [f.baseName.find(/(?<=chr)[A-Z0-9]+/), f] } // regex: chrX, chr11 -> X, 11 
      .set{ thousand_genomes }
    thousand_genomes_with_idx = tabix_chr(thousand_genomes)
    convert_imp5(thousand_genomes_with_idx)
    // joining requires same type, thousand genomes is a string (chrX)
    phased
      .map{ input -> tuple(input[0].toString(), input[1], input[2])}
      .set{phased_str_idx}
    impute_input = convert_imp5.out.chrom.join(phased_str_idx).combine(ch_map)
    impute5(impute_input)
    merge_imp(impute5.out.imputed.collect())
    // logs = fix_duplicates.out.log.concat(to_bcf.out.log).collect()
    // parse_logs("imputation", logs, "imputation_log.txt")

  emit:
    imputed_vcf = merge_imp.out.vcf

} 

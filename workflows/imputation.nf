// Imputation workflow
nextflow.preview.dsl = 2

// import modules
include {run_snpflip} from '../modules/popStrat.nf' // D1, reuse C4
include {flip_snps} from '../modules/popStrat.nf' // D2, D4, reuse C4
include {fix_duplicates} from '../modules/imputation.nf' // D3
include {to_bcf} from '../modules/imputation.nf' // D5 - D7
include {check_ref_allele} from '../modules/imputation.nf' // D8
include {bcf_to_vcf} from '../modules/imputation.nf' // D9 - D11
include {split_user_chrom} from '../modules/imputation.nf' // D12 - D13
include {phasing} from '../modules/imputation.nf' // D14 - D15
include {convert_imp5} from '../modules/imputation.nf' // D17 (D16 in container)
include {impute5} from '../modules/imputation.nf' // D18

// workflow component for snpqt pipeline
workflow imputation {
  take:
    ch_bed
    ch_bim
    ch_fam
    ch_db

  main:
    run_snpflip(ch_bed, ch_bim, ch_fam, ch_db)
    flip_snps(ch_bed, ch_bim, ch_fam, run_snpflip.out.rev, run_snpflip.out.ambig)
    fix_duplicates(flip_snps.out.bed, flip_snps.out.bim, flip_snps.out.fam)
    to_bcf(fix_duplicates.out.bed, fix_duplicates.out.bim, fix_duplicates.out.fam)
    check_ref_allele(to_bcf.out.bcf, ch_db)
    bcf_to_vcf(check_ref_allele.out.bcf)
    // OK, things start to get a bit complicated here
    // because we're working with individual chromosomes it's easiest to work with
    // tuples like [chr_num, chr.fa.gz, chr.fa.gz.csi, etc...]
    // that way we can join and combine files by chr_num as needed
    Channel.from(1..22).set{ chrom }
    split_user_chrom(bcf_to_vcf.out.vcf, bcf_to_vcf.out.idx, chrom)
    // maps for shapeit4
    Channel
      .fromPath("$baseDir/db/genetic_maps.b37.tar.gz", checkIfExists: true)
      .set { ch_map }
    user_chroms = split_user_chrom.out.chrom.combine(ch_map)
    phasing(user_chroms)
    // thousand genome reference data
    Channel
      .fromPath("$baseDir/db/*genotypes.vcf_updated.vcf.gz", checkIfExists: true)
      .map{ f -> [f.baseName.find(/\d+/).toInteger(), f] } // .toInteger() for join
      .set{ thousand_genomes }
    Channel
      .fromPath("$baseDir/db/*genotypes.vcf_updated.vcf.gz.csi", checkIfExists: true)
      .map{ f -> [f.baseName.find(/\d+/).toInteger(), f] } 
      .join( thousand_genomes )
      .set{ thousand_genomes_with_idx }
    convert_imp5(thousand_genomes_with_idx)
    impute5(convert_imp5.out.chrom.join(phasing.out.chrom).combine(ch_map))

  emit:
    imputed = impute5.out.imputed.collect()

} 
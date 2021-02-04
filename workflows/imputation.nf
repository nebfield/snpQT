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
include {phasing} from '../modules/imputation.nf' // D14 
include {bcftools_index_chr} from '../modules/imputation.nf' // D15
include {convert_imp5} from '../modules/imputation.nf' // D17 (D16 in container)
include {impute5} from '../modules/imputation.nf' // D18
include {parse_logs} from ../modules/qc.nf'

// workflow component for snpqt pipeline
workflow imputation {
  take:
    ch_bed
    ch_bim
    ch_fam

  main:
     Channel
      .fromPath("$baseDir/db/h37_squeezed.fasta", checkIfExists: true)
      .set { g37 }
    run_snpflip(ch_bed, ch_bim, ch_fam, g37)
    flip_snps(ch_bed, ch_bim, ch_fam, run_snpflip.out.rev, run_snpflip.out.ambig)
    fix_duplicates(flip_snps.out.bed, flip_snps.out.bim, flip_snps.out.fam)
    to_bcf(fix_duplicates.out.bed, fix_duplicates.out.bim, fix_duplicates.out.fam)
    Channel
      .fromPath("$baseDir/db/impute/All_20180423.vcf.gz", checkIfExists: true)
      .set{ dbsnp }
    Channel
      .fromPath("$baseDir/db/impute/All_20180423.vcf.gz.tbi", checkIfExists: true)
      .set{ dbsnp_idx }
    check_ref_allele(to_bcf.out.bcf, dbsnp, dbsnp_idx, g37)
    bcf_to_vcf(check_ref_allele.out.bcf)
    // OK, things start to get a bit complicated here
    // because we're working with individual chromosomes it's easiest to work with
    // tuples like [chr_num, chr.fa.gz, chr.fa.gz.csi, etc...]
    // that way we can join and combine files by chr_num as needed
    Channel.from(1..22).set{ chrom }
    split_user_chrom(bcf_to_vcf.out.vcf, bcf_to_vcf.out.idx, chrom)
    // maps for shapeit4
    Channel
      .fromPath("$baseDir/db/impute/genetic_maps.b37.tar.gz", checkIfExists: true)
      .set { ch_map }
    user_chroms = split_user_chrom.out.chrom.combine(ch_map)
    phasing(user_chroms)
    phased = bcftools_index_chr(phasing.out.chrom)
    // thousand genome reference data
    Channel
      .fromPath("$baseDir/db/impute/[0-9]*vcf.gz", checkIfExists: true)
      .map{ f -> [f.baseName.find(/\d+/).toInteger(), f] } // .toInteger() for join
      .set{ thousand_genomes }
    Channel
      .fromPath("$baseDir/db/impute/[0-9]*vcf.gz.csi", checkIfExists: true)
      .map{ f -> [f.baseName.find(/\d+/).toInteger(), f] } 
      .join( thousand_genomes )
      .set{ thousand_genomes_with_idx }
    convert_imp5(thousand_genomes_with_idx)
    impute5(convert_imp5.out.chrom.join(phased).combine(ch_map))
    logs = fix_duplicates.out.log.concat(to_bcf.out.log).collect()
    parse_logs("imputation", logs, "imputation.log")

  emit:
    imputed = impute5.out.imputed.collect()

} 
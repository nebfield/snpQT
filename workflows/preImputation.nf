// Imputation workflow
nextflow.enable.dsl = 2

// import modules
include {set_chrom_code} from '../modules/imputation.nf' // F1
include {run_snpflip} from '../modules/popStrat.nf' // F2, reuse D4
include {flip_snps} from '../modules/popStrat.nf' // F2, reuse D4
include {fix_duplicates} from '../modules/imputation.nf' // F3
include {to_vcf} from '../modules/imputation.nf' // F4
include {to_bcf} from '../modules/imputation.nf' // F5
include {check_ref_allele} from '../modules/imputation.nf' // F6
include {bcf_to_vcf} from '../modules/imputation.nf' // F7
include {parse_logs} from '../modules/qc.nf' // E13

// workflow component for snpqt pipeline
workflow preImputation {
  take:
    ch_bed
    ch_bim
    ch_fam

  main:
    set_chrom_code(ch_bed, ch_bim, ch_fam)
    Channel
      .fromPath("${params.db}/h37_squeezed.fasta", checkIfExists: true)
      .set { g37 }
    run_snpflip(set_chrom_code.out.bed, set_chrom_code.out.bim, set_chrom_code.out.fam, g37)
    flip_snps(ch_bed, ch_bim, ch_fam, run_snpflip.out.rev, run_snpflip.out.ambig)
    fix_duplicates(flip_snps.out.bed, flip_snps.out.bim, flip_snps.out.fam)
    to_vcf(fix_duplicates.out.bed, fix_duplicates.out.bim, fix_duplicates.out.fam)
	to_bcf(to_vcf.out.vcf)
    Channel
      .fromPath("${params.db}/All_20180423.vcf.gz", checkIfExists: true)
      .set{ dbsnp }
    Channel
      .fromPath("${params.db}/All_20180423.vcf.gz.tbi", checkIfExists: true)
      .set{ dbsnp_idx }
    check_ref_allele(to_bcf.out.bcf, dbsnp, dbsnp_idx, g37)
    bcf_to_vcf(check_ref_allele.out.bcf)
    // logs = fix_duplicates.out.log.concat(to_vcf.out.log).collect()
    // parse_logs("pre_imputation", logs, "pre_imputation_log.txt")

  emit:
    vcf = bcf_to_vcf.out.vcf
}  

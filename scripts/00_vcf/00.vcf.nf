/* VCF sanity checking nextflow pipeline
   Parameters:
     params.infile:
     params.outdir:
     params.ref_fasta 
*/

log.info """\
         snpQT: VCF sanity checking 
         input file: ${params.infile}
         output directory: ${params.outdir}
         reference fasta: ${params.ref_fasta}
         """
         .stripIndent()

Channel
    .fromPath( params.infile )
    .ifEmpty { error "Cannot find: ${params.infile}" }
    .set { in_file } 

// STEP A1: map to hg37 (TODO) ------------------------------------------------

// STEP A2: Convert VCF to PLINK format ---------------------------------------
process convert_VCF {
    echo true
    container 'snpqt'
    publishDir "$baseDir/results", mode: 'copy', overwrite: true, 
      pattern: "*.pdf"

    input:
    file in_file

    output:
    file "dataset_2.bim" into dataset_2_bim
    file "dataset_2.bed" into dataset_2_bed
    file "dataset_2.fam" into dataset_2_fam
    file "dataset_2*" into dataset_2

    """
    plink --vcf $in_file --make-bed --keep-allele-order --out dataset_2 
    """
}

// STEP A3: Flip SNPs with snpflip --------------------------------------------
// TODO: move this to pre-imputation QC
reference_fasta = Channel.fromPath(params.ref_fasta).collect()
    
process flip_snp {
    echo true
    container 'snpqt'

    input:
    file dataset_2_bim
    file reference_fasta

    output:
    file "snpflip.reverse" into snpflip_reverse
    file "snpflip.ambiguous"
    file "snpflip.annotated_bim"

    """
    snpflip -b $dataset_2_bim -f $reference_fasta -o snpflip
    """
}

process SNP_flip_PLINK {
    echo true
    container 'snpqt'

    input:
    file dataset_2 
    file snpflip_reverse

    output:
    file "dataset_3*" into dataset_3

    """
    plink --bfile dataset_2 --flip $snpflip_reverse --make-bed --out dataset_3
    """
 }

// STEP A4: Remove ambiguous SNPs (TODO) --------------------------------------

// STEP A5: Remove duplicated SNPs (TODO) -------------------------------------

// STEP A6: Convert back to VCF (TODO: combine, plink can write bgzip)
process plink_to_vcf {
    echo true
    container 'snpqt'

    input:
    file dataset_3

    output:
    file "dataset_4.vcf" into dataset_4_vcf
    file "dataset_4.log" 
    file "dataset_4.nosex"

    """
    plink --bfile dataset_3 --recode vcf --keep-allele-order --out dataset_4
    """
}

process vcf_gz {
    echo true
    container 'snpqt'

    input:
    file dataset_4_vcf

    output:
    file "dataset_4.vcf.gz" into dataset_5

    """
    bgzip -c $dataset_4_vcf > dataset_4.vcf.gz
    """
}

process vcf_index {
    echo true
    container 'snpqt'

    input:
    file dataset_5

    """
    bcftools index $dataset_5
    """
}
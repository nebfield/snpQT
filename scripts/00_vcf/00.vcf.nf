/* VCF sanity checking nextflow pipeline
   Parameters:
     params.infile:
     params.outdir:
     params.ref_fasta 
*/

log.info """\
         snpQT: VCF sanity checking 
         input file: ${params.infile}
         """
         .stripIndent()

Channel
    .fromPath( params.infile )
    .ifEmpty { error "Cannot find: ${params.infile}" }
    .set { in_file } 

// STEP A1: Change the chr ids (TODO) ------------------------------------------

// STEP A2: Run liftOver to map genome build -----------------------------------

// STEP A3: Change chr ids (TODO) ----------------------------------------------

// STEP A4: Convert VCF to PLINK format ----------------------------------------
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
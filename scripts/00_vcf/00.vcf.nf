/* VCF sanity checking nextflow pipeline
   Parameters:
     params.infile:
     params.outdir:
     params.ref_fasta 
*/

log.info """\
         snpQT: VCF sanity checking 
         """
         .stripIndent()

// STEP A1: Change the chr ids (TODO) ------------------------------------------


// STEP A2: Run liftOver to map genome build -----------------------------------

// STEP A3: Change chr ids (TODO) ----------------------------------------------

// STEP A4: Convert VCF to PLINK format ----------------------------------------

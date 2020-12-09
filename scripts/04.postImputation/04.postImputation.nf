// Post-imputation (IN PROGRESS)
// ============================================================================

params.inimp = "../../results/imputation/imputed_chr*"
params.outdir = "$baseDir/../../results"
outdir = params.outdir + '/postImputation/'

log.info """\
         snpQT step04: Post-imputation QC 
         """
         .stripIndent()

// User's imputed chromosomes --------------------------------------------------
// extract chromosome number from file name to make a tuple of [chr, file path]
Channel
  .fromPath( params.inimp )
  .first() // for testing (TODO: REMOVE) =======================================
  .map{ f -> [f.baseName.find(/\d+/), f] }
  .set { inimp }

// Pre-imputation 
// =============================================================================

// STEP E1: Filter all poorly imputed variants based on info score
// (check impute5 output)
process filter {
    container 'snpqt'

    input:
    tuple val(chr_no), file(chr) from inimp

    output:
    tuple val(chr_no), file("E1.vcf.gz") into E1
    
    shell:
    '''
    bcftools view -Ou -i 'INFO>0.7' -q 0.05:minor !{chr} -Oz -o E1.vcf.gz
    '''
}

// STEP E2: Annotate missing variant ids and convert .vcf.gz to binary plink
process annotate {
    container 'snpqt'

  input:
  tuple val(chr_no), file(chr) from E1

  output:
  file "E2.*" into E2
  
  shell:
  '''
  plink2 --vcf !{chr} \
    --vcf-idspace-to _ \
    --const-fid \
    --set-missing-var-ids @:#:\\$r:\\$a \
    --new-id-max-allele-len 1000 \
    --make-bed \
    --out E2
  '''
}
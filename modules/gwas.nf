// Genomw Wide Association Study (GWAS) module
// ============================================================================
// Input: B11 genomic output and C10 output (covar_pca)

params.inbed = "../../results/sample_qc/sample_variant_qc.*"
params.inbim = "../../results/sample_qc/sample_variant_qc.bim"
params.infam = "../../results/sample_qc/sample_variant_qc.fam"
params.outdir = "$baseDir/../../results"
outdir = params.outdir + '/gwas/'

log.info """\
         snpQT step03: Imputation
         input file: ${params.inbim}
         outdir: ${params.outdir}
         """
         .stripIndent()

// User's plink files ---------------------------------------------------------
Channel.fromPath( params.inbed ).into { inbed } 
Channel.fromPath( params.inbim ).into { inbim } 
Channel.fromPath( params.infam ).into { infam } 

// Step F1: Run logistic regression, adjusting for covariates

process gwas {
    publishDir outdir, mode: 'copy', overwrite: true
    
    input:
    file inbed
    file inbim
    file infam

    output:
    file 'logistic_results*'
    
    shell:
    '''
    plink --bfile !{inbed.baseName} \
        --covar covar_pca.txt \
	--ci 0.95 \
	--logistic \
	--allow-no-sex \
	--out logistic_results
    '''
}
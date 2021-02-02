// Step F1: Run logistic regression, adjusting for covariates

process run_gwas {    
    input:
    path(bed)
    path(bim)
    path(fam)
    path(covar)

    output:
    path "*.logistic", emit: logistic
    
    shell:
    '''
    plink --bfile !{bed.baseName} \
      --covar !{covar} \
      --ci 0.95 \
      --logistic \
      --allow-no-sex \
      --out logistic_results
    plink --bfile !{bed.baseName} \
      --ci 0.95 \
      --logistic \
      --allow-no-sex \
      --out logistic_results_nocovars
    '''
}

process plot {
    publishDir "${params.results}/gwas/$id/", mode: 'copy'
    
    input:
    tuple id, path(logistic)

    output:
    path "qqplot.png"
    
    shell:
    '''
    qqplot.R !{logistic}
    manhattan.R !{logistic}
    '''
}

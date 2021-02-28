// Step F1: Run logistic regression, adjusting for covariates

process run_gwas {    
    publishDir "${params.results}/gwas/files/", mode: 'copy'

    input:
    path(bed)
    path(bim)
    path(fam)
    path(covar)

    output:
    path "*.logistic", emit: logistic
    path "*.log", emit: log
    
    shell:
    '''
    plink --bfile !{bed.baseName} \
      --covar !{covar} \
      --ci 0.95 \
      --logistic hide-covar\
      --allow-no-sex \
      --out logistic_results
    plink --bfile !{bed.baseName} \
      --ci 0.95 \
      --logistic hide-covar\
      --allow-no-sex \
      --out logistic_results_nocovars
    '''
}

process plot {
    publishDir "${params.results}/gwas/figures/", mode: 'copy'
    
    input:
    tuple id, path(logistic)

    output:
    path "*_qqplot.png", emit: qqplot
    path "*_manhattan.png", emit: manhattan
 
    shell:
    '''
    qqplot.R !{logistic}
    manhattan.R !{logistic}
    cp qqplot.png !{id}_qqplot.png
    cp manhattan.png !{id}_manhattan.png
    '''
}

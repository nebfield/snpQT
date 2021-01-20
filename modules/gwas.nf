// Step F1: Run logistic regression, adjusting for covariates

process run_gwas {  
    input:
    path(bed)
    path(bim)
    path(fam)
    path(covar)

    output:
    path "logistic_results.assoc.logistic", emit: logistic
    
    shell:
    '''
    plink --bfile !{bed.baseName} \
        --covar !{covar} \
	--ci 0.95 \
	--logistic \
	--allow-no-sex \
	--out logistic_results
    '''
}

process plot {
    input:
    path(logistic)

    output:
    path "qqplot.png", emit: qqplot
    
    shell:
    '''
    qqplot.R !{logistic}
    # manhattan.R !{logistic} # qqman broken >:( 
    '''
}
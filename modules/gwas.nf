// Step I1: Run logistic regression, adjusting for covariates
process run_gwas {    
    publishDir "${params.results}/gwas/files/", mode: 'copy'

    input:
    path(bed)
    path(bim)
    path(fam)
    path(covar)

    output:
    path "*.logistic.hybrid", emit: logistic
    path "*.log", emit: log
    
    shell:
    '''
    plink2 --bfile !{bed.baseName} \
      --covar !{covar} \
      --ci 0.95 \
      --glm hide-covar firth-fallback\
	  --output-chr 26 \
      --out logistic_results
    plink2 --bfile !{bed.baseName} \
      --ci 0.95 \
      --glm hide-covar firth-fallback\
	  --output-chr 26 \
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
    # NA p values sometimes happen
    awk '!/'NA'/' !{logistic} > no_na.logistic
	# Remove hash from beginning of line
	sed 's/#//' no_na.logistic | sed 's/LOG(OR)_SE/SE/' | sed 's/CHROM/CHR/'|  sed 's/POS/BP/' | sed 's/ID/SNP/' > valid.logistic
	# Run plots
    qqplot.R valid.logistic
    manhattan.R valid.logistic
    cp qqplot.png !{id}_qqplot.png
    cp manhattan.png !{id}_manhattan.png
    '''
}

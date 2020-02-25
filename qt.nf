
params.infile = "../data/als.vcf"
params.outdir = "$baseDir/results"
params.famfile = "../data/plinkForDbgap12319.fam"
params.sex_impute = false
params.highldregion = "../data/highldregion_37.txt" // TODO: update me 
params.ref_file = "$baseDir/../data/1kG_MDS3*"

log.info """\
         snpQT: make your SNPs cute 
         input file: ${params.infile}
         outdir: ${params.outdir}
         fam file: ${params.famfile}
         sex impute: ${params.sex_impute}
         """
         .stripIndent()

Channel
    .fromPath( params.infile )
    .ifEmpty { error "Cannot find: ${params.infile}" }
    .set { in_file } 

Channel
    .fromPath( params.famfile )
    .ifEmpty { error "Cannot find: ${params.famfile}" }
    .set { fam_file } 

// vcf to bed + fam https://www.biostars.org/p/313943/
process check_sex {
    echo true
    container 'snpqt'
    publishDir "$baseDir/results", mode: 'copy', overwrite: true, 
      pattern: "*.pdf"

    input:
    file in_file
    file fam_file

    output:
    file "plink.sexcheck" into sexcheck 
    file "data_clean.bed" into sex_checked_bed
    file "data_clean.bim" into sex_checked_bim
    file "data_clean.fam" into sex_checked_fam

    """
    plink --make-bed --vcf $in_file  --out data &>/dev/null
    # the input stage will change as the pipeline is developed 
    echo 'Pipeline input: ' && grep 'pass' data.log
    cp $fam_file data.fam
    plink --bfile data --check-sex &>/dev/null

    # Remove samples with ambiguous sex phenotypes
    if [ -f plink.nosex ]; then
        plink --bfile data --remove plink.nosex --make-bed \
          --out data_ambig &>/dev/null
        mv data_ambig.bed data.bed
        mv data_ambig.bim data.bim
        mv data_ambig.fam data.fam
    fi

    # Identify the samples with sex discrepancy 
    grep "PROBLEM" plink.sexcheck | awk '{print \$1,\$2}'> \
      problematic_samples.txt
    # Delete all problematic samples
    plink --bfile data --remove problematic_samples.txt --make-bed \
      --out data_clean &>/dev/null
    echo 'Check sex output: ' && grep 'pass' data_clean.log
    """
}

process sex_impute {
    echo true

    when:
    params.sex_impute == true

    """
    echo 'sex impute true!'
    """ 
}

process plot_sex {
    publishDir "$baseDir/results", mode: 'copy', overwrite: true, 
      pattern: "*.png"
    container 'rocker/tidyverse:3.6.1' 

    input:
    file sexcheck 

    output:
    file "*.pdf"

    """
    gender_check.R 
    """
}

process extract_autosomal {
    echo true
    container 'snpqt'

    input:
    file sex_checked_bed   
    file sex_checked_bim
    file sex_checked_fam

    output:
    file "autosomal.*" into autosomal

    """
    # Extract only autosomal chromosomes (for studies that don't want to 
    # include sex chromosomes)
    awk '{ if (\$1 >= 1 && \$1 <= 22) print \$2 }' $sex_checked_bim > \
      autosomal_SNPs.txt 
    plink --bfile data_clean --extract autosomal_SNPs.txt --make-bed \
      --out autosomal &>/dev/null
    echo 'Extract autosomal output:' && grep 'pass' autosomal.log
    """ 
}

process missing {
    echo true
    container 'snpqt'

    input:
    file autosomal 
    
    output:
    file "plink.imiss" into imiss_before, imiss_relatedness
    file "plink.lmiss" into lmiss_before
    file "cleaned*" into missing, missing_het
 
    """
    plink --bfile autosomal --missing --out plink &>/dev/null
    plink --bfile autosomal --geno 0.1 --mind 0.1 --make-bed \
      --out cleaned &>/dev/null
    echo 'Check missing output:' && grep 'pass' cleaned.log 
    """
}

process plot_missing {
    publishDir "$baseDir/results", mode: 'copy', overwrite: true, 
      pattern: "*.png"
    container 'rocker/tidyverse:3.6.1' 
 
    input:
    file imiss_before
    file lmiss_before
 
    output:
    file "sample_callrate.png"
    file "variant_callrate.png"
 
    """
    plot_missing.R $imiss_before $lmiss_before
    """
} 

Channel
    .fromPath( params.highldregion )
    .ifEmpty { error "Cannot find: ${params.highldregion}" }
    .set { high_ld_file } 

process heterozygosity_rate {
    container 'snpqt'

    input:
    file high_ld_file 
    file missing

    output:
    file "only_indep_snps*" into plot_het
    file "independent_SNPs.prune.in" into ind_SNPs, ind_SNPs_popstrat

    """
    plink --bfile cleaned --exclude $high_ld_file --indep-pairwise 50 5 0.2 \
      --out independent_SNPs --range 
    plink --bfile cleaned --extract independent_SNPs.prune.in --het \
      --out only_indep_snps 
    """
}

process plot_heterozygosity { 
    publishDir "$baseDir/results", mode: 'copy', overwrite: true, 
      pattern: "*.png"
    container 'rocker/tidyverse:3.6.1' 

    input: 
    file plot_het

    output:
    file "het_failed_samples.txt" into het_failed
    file "heterozygosity_rate.png"

    """
    plot_heterozygosity.R
    """
}

process heterozygosity_prune {
    echo true
    container 'snpqt'

    input:
    file missing_het
    file het_failed

    output:
    file "het_pruned.*" into het_pruned

    """
    cut -f 1,2 $het_failed > het_failed_plink.txt
    plink --bfile cleaned --make-bed --out het_pruned --remove \
      het_failed_plink.txt &>/dev/null
    echo 'Check heterozygosity rate:' && grep 'pass' het_pruned.log
    """
}

process relatedness {
    echo true
    container 'snpqt'
  
    input:
    file het_pruned
    file ind_SNPs
    file imiss_relatedness

    output:
    file "pihat_pruned*" into relatedness

    """
    plink --bfile het_pruned --extract $ind_SNPs --genome --min 0.125 \
      --out pihat_0.125 &>/dev/null
    # Identify all pairs of relatives with pihat > 0.125 and exclude one of the
    # relatives of each pair, having the most missingness. Output those failing
    # samples to pihat_failed_samples.txt
    run_IBD_QC.pl &>/dev/null
    plink --bfile het_pruned --remove pihat_failed_samples.txt \
      --make-bed --out pihat_pruned &>/dev/null
    echo 'Check relatedness:' && grep 'pass' pihat_pruned.log
    """
}

process missing_phenotype {
    echo true
    container 'snpqt'

    input:
    file relatedness

    output:
    file "missing*" into missing_pheno, missing_pop_strat

    """
    plink --bfile pihat_pruned --prune --make-bed --out missing &>/dev/null
    echo 'Remove missing phenotypes:' && grep 'pass' missing.log 
    """
}

thousand_ref_file = Channel.fromPath(params.ref_file).collect()

process download_pop_info {
    echo true
    container 'snpqt'
    
    output:
    file "race_1kG.txt" into racefile_1kG

    """
    download_racefile.sh
    """
}

process prep_pop_strat {
    echo true
    container 'snpqt'
    stageInMode 'copy'

    input:
    file thousand_ref_file
    file missing_pheno
    file racefile_1kG

    output:
    file "MDS_merge*" into pop_strat_mds
    file "racefile.txt" into racefile 

    """
    prep_pop_strat.sh 
    """
}

process mds_pop_strat {
    echo true
    container 'snpqt'

    input:
    file pop_strat_mds
    file ind_SNPs_popstrat
    file missing_pop_strat

    output:
    file "MDS_merge*" into mds
    file "plink_6*" into hardy_in
    file "covar_mds.txt" into covar_mds

    """
    mds.sh
    """
}

process plot_mds {
    echo true
    container 'rocker/tidyverse:3.6.1' 
    publishDir "$baseDir/results", mode: 'copy', overwrite: true, 
        pattern: "*.png"

    input:
    file racefile
    file mds

    output:
    file "MDS.png"

    """
    plot_MDS.R 
    """
}

process hardy {
    echo true
    container 'snpqt'

    input:
    file hardy_in

    output:
    file "*.hwe" into hwe
    file "plink_7*" into hardy_out

    """
    hardy.sh &>/dev/null
    grep 'pass' plink_7.log
    """
}

process plot_hardy {
    echo true
    container 'rocker/tidyverse:3.6.1' 
    publishDir "$baseDir/results", mode: 'copy', overwrite: true, 
        pattern: "*.png"

    input:
    file hwe
    
    output:
    file "*.png" 

    """
    hwe.R
    """
}

process maf {
    echo true
    container 'snpqt'

    input:
    file hardy_out

    output:
    file "plink_8*" into maf_output

    """
    plink --bfile plink_7 --freq --out MAF_check &>/dev/null

    # Remove SNPs with a low MAF frequency
    plink --bfile plink_7 --maf 0.05 --make-bed --out plink_8 &>/dev/null
    grep 'pass' plink_8.log
    """
}

process plot_maf {
    echo true
    container 'rocker/tidyverse:3.6.1' 

    """
    # TODO: missing MAF_check.R
    """
}

process case_control_status {
    echo true
    container 'snpqt'

    input:
    file maf_output

    output:
    file "plink_9*" into case_control_status_output

    """
    case_control_status.sh &>/dev/null
    grep 'pass' plink_9.log
    """
}

process log_reg {
    echo true
    container 'snpqt'

    input:
    file case_control_status_output
    file covar_mds

    """
    ## adjusting for covariates doesn't work :/ 

    plink --bfile plink_9 --covar covar_mds.txt --logistic --ci 0.95 \
        --out logistic_results

    plink --bfile plink_9 --logistic --ci 0.95 --out logistic_results1

    plink --bfile plink_9 --covar covar_mds.txt --logistic --hide-covar \
        --out logistic_results2

    # Note, we use the option -–hide-covar to only show the additive 
    # results of the SNPs in the output file.
    # Remove NA values, those might give problems generating plots in later
    # steps.

    awk '!/'NA'/' logistic_results.assoc.logistic \
        > logistic_results.assoc_3.logistic

    # TODO: missing qqPlot.R
    # R --vanilla --slave --args logistic_results.assoc.logistic ADD 
    # qq_1.pdf< qqPlot.R
    """
}

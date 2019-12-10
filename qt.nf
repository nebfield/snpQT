
params.infile = "$baseDir/data/wp2_good.vcf.gz"
params.outdir = "$baseDir/results"
params.famfile = "$baseDir/data/wp2_good.fam"
params.sex_impute = false
params.highldregion = "$baseDir/data/highldregion_37.txt" // TODO: update me 

log.info """\
         snpQT: make your SNPs cute 
         input file : ${params.infile}
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
    container 'quay.io/biocontainers/plink:1.90b6.12--heea4ae3_0'
    publishDir "$baseDir/results", mode: 'copy', overwrite: true, 
      pattern: "*.pdf"

    input:
    file in_file
    file fam_file

    output:
    file "*.pdf"
    file "data_clean.bed" into sex_checked_bed
    file "data_clean.bim" into sex_checked_bim
    file "data_clean.fam" into sex_checked_fam

    """
    plink --make-bed --vcf $in_file  --out data &>/dev/null
    # the input stage will change as the pipeline is developed 
    echo 'Pipeline input: ' && grep 'pass' data.log
    cp $fam_file data.fam
    plink --bfile data --check-sex &>/dev/null

    gender_check.R &>/dev/null

    # Identify the samples with sex discrepancy 
    grep "PROBLEM" plink.sexcheck | awk '{print \$1,\$2}'> \
      problematic_samples.txt
    
    # Delete all problematic samples
    plink --bfile data --remove problematic_samples.txt --make-bed \
      --out data_clean &>/dev/null
    echo 'Check sex output: ' && grep 'pass' data_clean.log

    # Remove samples with ambiguous sex phenotypes
    if [ -f plink.nosex ]; then
        plink --bfile data_clean --remove plink.nosex --make-bed \
          --out data_ambig
        mv data_ambig.bed data_clean.bed
        mv data_ambig.bim data_clean.bim
        mv data_ambig.fam data_clean.fam
        echo 'Check sex output:' && grep 'pass' data_ambig.log
    fi
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

process extract_autosomal {
    echo true
    container 'quay.io/biocontainers/plink:1.90b6.12--heea4ae3_0'

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
    container 'quay.io/biocontainers/plink:1.90b6.12--heea4ae3_0'

    input:
    file autosomal 
    
    output:
    file "plink.imiss" into imiss_before
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
    container 'quay.io/biocontainers/plink:1.90b6.12--heea4ae3_0'

    input:
    file high_ld_file 
    file missing

    output:
    file "only_indep_snps*" into plot_het
    file "independent_SNPs.prune.in" into ind_SNPs

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
    container 'quay.io/biocontainers/plink:1.90b6.12--heea4ae3_0'

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
    container 'quay.io/biocontainers/plink:1.90b6.12--heea4ae3_0'
  
    input:
    file het_pruned
    file ind_SNPs

    """
    plink --bfile het_pruned --extract $ind_SNPs --genome --min 0.125 \
      --out pihat_0.125 &>/dev/null
    """
}

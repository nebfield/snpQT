// TODO: make sure parameters exist
// TODO: default sensible parameters? maybe a test data directory?
// TODO: single log file

// params.inbed = "../data/als_sub.vcf.gz"
// params.inbim = "../data/als_sub.vcf.gz"
// params.infam = "../data/subset.fam"
params.outdir = "$baseDir/../../results"
outdir = params.outdir + '/sample_qc/'

log.info """\
         snpQT step01: sample variant quality control
         input file: ${params.inbed}
         outdir: ${params.outdir}
         fam file: ${params.infam}
         """
         .stripIndent()

Channel.fromPath( params.inbed ).set { in_bed }
Channel.fromPath( params.inbim ).set { in_bim } 
Channel.fromPath( params.infam ).set { in_fam } 

Channel.fromPath("$SNPQT_DB_DIR/PCA.exclude.regions.b37.txt").set { exclude_regions } 


// STEP B1: Remove SNPs < 90% missingness --------------------------------------
process missingness {
  input:
  file in_bed
  file in_bim
  file in_fam

  output:
  file "missingness.log" into missingness_logs
  file "plink_1*" into missingness_bfiles

  shell:
  '''  
  cp !{in_fam} dataset_4.fam # rename fam file to match input bim / bed
  plink --make-bed --bfile dataset_4 --out data &>/dev/null 
  
  echo 'Pipeline input: ' && grep 'pass' data.log > log.txt
  cp !{in_fam} data.fam
  echo 'Step B1: missingness' >> log.txt
  plink --bfile data --geno 0.1 --make-bed  --out plink_1 &>/dev/null
  echo 'Missingness output: ' && grep 'pass' plink_1.log
  cat *.log > missingness.log
  '''
}

// STEP B2: Check missingness rate ---------------------------------------------
// TODO: user set threshold
process plot_missingness {
  echo true
  publishDir outdir, mode: 'copy', overwrite: true, pattern: "*.png"

  input:
  file missingness_bfiles

  output:
  file "sample_missingness.png"
  file "plink.imiss" into imiss_relatedness
  file "plink_3*" into missingness_bfiles_pruned

  shell:
  '''
  plink --bfile plink_1 \
    --missing \
    --out plink_2 &>/dev/null
  plot_sample_missingness.R plink_2.imiss
  # TODO user set mind threshold 
  plink --bfile plink_1  \
    --make-bed \
    --mind 0.02 \
    --out plink_3 &>/dev/null
  mv plink_2.imiss plink.imiss # run_IBD_QC.pl
  '''
}

// STEP B3: Remove samples with sex mismatch -----------------------------------
process check_sex {
    input:
    file missingness_bfiles_pruned

    output:
    file "plink.sexcheck" into sexcheck 
    file "data_clean.bed" into sex_checked_bed
    file "data_clean.bim" into sex_checked_bim
    file "data_clean.fam" into sex_checked_fam

    // vcf to bed + fam https://www.biostars.org/p/313943/
    shell:
    '''
    plink --bfile plink_3 \
      --check-sex 
    # Identify the samples with sex discrepancy 
    grep "PROBLEM" plink.sexcheck | awk '{print $1,$2}'> \
      problematic_samples.txt
    # Delete all problematic samples
    plink --bfile plink_3 \
      --remove problematic_samples.txt \
      --make-bed \
      --out data_clean 
    echo 'Check sex output: ' && grep 'pass' data_clean.log
    '''
}

process plot_sex {
    publishDir outdir, mode: 'copy', overwrite: true

    input:
    file sexcheck 

    output:
    file "sexcheck.png"

    shell:
    '''
    plot_sex.R !{sexcheck}
    '''
}

// STEP B4: Remove sex chromosomes ---------------------------------------------
// shell for awk ($ confuses nextflow)
process extract_autosomal {
    echo true

    input:
    file sex_checked_bed   
    file sex_checked_bim
    file sex_checked_fam

    output:
    file "autosomal.*" into autosomal, autosomal_het

    shell:
    '''
    # Extract only autosomal chromosomes (for studies that don't want to 
    # include sex chromosomes)
    awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' !{sex_checked_bim} > \
      autosomal_SNPs.txt 
    plink --bfile data_clean \
      --extract autosomal_SNPs.txt \
      --make-bed \
      --out autosomal &>/dev/null
    echo 'Extract autosomal output:' && grep 'pass' autosomal.log
    ''' 
}

// STEP B5: Remove SNPs with extreme heterozygosity ----------------------------
// Extract highly independent SNPs based on LD and remove MHC region 

process heterozygosity_rate {
    input:
    file exclude_regions 
    file autosomal

    output:
    file "only_indep_snps*" into plot_het
    file "independent_SNPs.prune.in" into ind_SNPs, ind_SNPs_popstrat

    shell:
    '''
    plink --bfile autosomal \
      --exclude 1{exclude_regions} \
      --indep-pairwise 50 5 0.2 \
      --out independent_SNPs \
      --range
    plink --bfile autosomal \
      --extract independent_SNPs.prune.in \
      --het \
      --out only_indep_snps 
    '''
}

process plot_heterozygosity { 
    publishDir outdir, mode: 'copy', overwrite: true, pattern: "*.png"

    input: 
    file plot_het

    output:
    file "het_failed_samples.txt" into het_failed
    file "heterozygosity_rate.png"

    shell:
    '''
    # plot and get outliers list 
    plot_heterozygosity.R
    '''
}

process heterozygosity_prune {
    echo true

    input:
    file autosomal_het
    file het_failed

    output:
    file "het_pruned.*" into het_pruned

    shell:
    '''
    cut -f 1,2 !{het_failed} > het_failed_plink.txt
    plink --bfile autosomal \
      --make-bed \
      --out het_pruned \
      --remove \
      het_failed_plink.txt &>/dev/null
    echo 'Check heterozygosity rate:' && grep 'pass' het_pruned.log
    '''
}

// STEP B6: Remove relatives ---------------------------------------------------
process relatedness {
    echo true
    container 'snpqt'
  
    input:
    file het_pruned
    file ind_SNPs
    file imiss_relatedness

    output:
    file "pihat_pruned*" into relatedness

    shell:
    '''
    plink --bfile het_pruned \
      --extract $ind_SNPs \
      --genome --min 0.125 \
      --out pihat_0.125 &>/dev/null
    # Identify all pairs of relatives with pihat > 0.125 and exclude one of the
    # relatives of each pair, having the most missingness. Output those failing
    # samples to pihat_failed_samples.txt
    run_IBD_QC.pl &>/dev/null
    plink --bfile het_pruned \
      --remove pihat_failed_samples.txt \
      --make-bed \
      --out pihat_pruned &>/dev/null
    echo 'Check relatedness:' && grep 'pass' pihat_pruned.log
    '''
}

// STEP B7: Remove samples with missing phenotypes -----------------------------
process missing_phenotype {
    echo true
    container 'snpqt'

    input:
    file relatedness

    output:
    file "missing*" into missing_pheno

    shell:
    '''
    plink --bfile pihat_pruned \
      --prune \
      --make-bed \
      --out missing &>/dev/null
    echo 'Remove missing phenotypes:' && grep 'pass' missing.log 
    '''
}

// STEP B8: Missingness per variant --------------------------------------------
// TODO: user defined?

process missingness_per_variant {
    publishDir outdir, mode: 'copy', overwrite: true, pattern: "*.png"

    input:
    file missing_pheno

    output:
    file "*.png"
    file "mpv.*" into mpv

    shell:
    '''
    plink --bfile missing --missing 
    plot_variant_missingness.R plink.lmiss
    plink --bfile missing --geno 0.05 --make-bed --out mpv 
    '''
}

// STEP B9: Hardy_Weinberg equilibrium (HWE) -----------------------------------
// shell for awk ($ confuses nextflow)
process hardy {
  publishDir outdir, mode: 'copy', overwrite: true, pattern: "*.png"

  input:
  file mpv

  output:
  file "hwe.*" into hwe
  file "*.png"

  shell: 
  '''
  plink --bfile mpv --hardy &>/dev/null
  # sample 1% of SNPs
  head -n1 plink.hwe > plink_sub.hwe
  perl -ne 'print if (rand() < 0.01)' <(tail -n +2 plink.hwe) >> plink_sub.hwe
  awk '{ if ($3=="TEST" || $3=="UNAFF" && $9 <0.001) print $0 }' \
	  plink.hwe > plinkzoomhwe.hwe
  hwe.R plink_sub.hwe ""
  hwe.R plinkzoomhwe.hwe "strongly deviating SNPs only"
  plink --bfile mpv --hwe 1e-7 --make-bed --out hwe &>/dev/null
  '''
}

// STEP B10: Remove low minor allele frequency (MAF) ---------------------------
process maf {
  publishDir outdir, mode: 'copy', overwrite: true, pattern: "*.png"

  input:
  file hwe

  output:
  file "MAF_check.*" into maf_check
  file "maf.png"

  shell:
  '''
  plink --bfile hwe \
    --freq \
    --out MAF_check
  plot_maf.R MAF_check.frq
  plink --bfile hwe \
    --maf 0.05 \
    --make-bed \
    --out MAF_check
  '''
}

// STEP B11: Test missingness in case / control status -------------------------
// shell for awk ($ confuses nextflow)
process test_missing {
  publishDir outdir, mode: 'copy', overwrite: true, pattern: "*.bed"
  publishDir outdir, mode: 'copy', overwrite: true, pattern: "*.bim"
  publishDir outdir, mode: 'copy', overwrite: true, pattern: "*.fam"

  input:
  file maf_check

  output:
  file "sample_variant_qc*"

  shell:
  '''
  plink --bfile MAF_check \
    --test-missing
  awk '{ if ($5 < 10e-7) print $2 }' plink.missing > fail_missingness.txt
  plink --bfile MAF_check \
    --exclude fail_missingness.txt \
    --make-bed \
    --out sample_variant_qc
  '''
}

// Finished! ------------------------------------------------------------------

// qc processes

// STEP B1: Remove SNPs < 90% missingness --------------------------------------
process variant_missingness {
  input:
  path(in_bed)
  path(in_bim)
  path(in_fam)

  output:
  path "B1.bed", emit: bed
  path "B1.bim", emit: bim
  path "B1.fam", emit: fam

  shell:
  '''  
  cp !{in_fam} !{in_bed.baseName}.fam # rename fam file to match input bim / bed
  plink --bfile !{in_bed.baseName} \
      --make-bed \
      --out data 
  plink --bfile data \
      --geno 0.1 \
      --make-bed  \
      --out B1 
  '''
}

// STEP B2: Check missingness rate ---------------------------------------------
process individual_missingness {
  input:
  path(B1_bed)
  path(B1_bim)
  path(B1_fam)

  output:
  path "B2.bed", emit: bed
  path "B2.bim", emit: bim
  path "B2.fam", emit: fam
  path "missing.imiss", emit: imiss
  
  shell:
  '''
  plink --bfile !{B1_bed.baseName} \
    --missing \
    --out missing 
  
  plink --bfile !{B1_bed.baseName} \
    --make-bed \
    --mind !{params.mind} \
    --out B2 
  '''
}

process plot_missingness {
  publishDir "${params.results}/qc/", mode: 'copy'

  input:
  path(missing_imiss)

  output:
  path "sample_missingness.png", emit: figure

  shell:
  '''
  plot_sample_missingness.R !{missing_imiss}
  '''
}

// STEP B3: Remove samples with sex mismatch -----------------------------------
// --check-sex requires at least one X chromosome so it has to be completed 
// before excluding non-automosomal SNPs 

process check_sex {
    input:
    path(B2_bed)
    path(B2_bim)
    path(B2_fam)
    
    output:
    path "B3.bed", emit: bed 
    path "B3.bim", emit: bim
    path "B3.fam", emit: fam
    path "plink.sexcheck", emit: sexcheck
    
    // vcf to bed + fam https://www.biostars.org/p/313943/
    shell:
    '''
    plink --bfile !{B2_bed.baseName} \
      --check-sex 
    # Identify the samples with sex discrepancy 
    grep "PROBLEM" plink.sexcheck | awk '{print $1,$2}'> \
      problematic_samples.txt
    # Delete all problematic samples
    plink --bfile !{B2_bed.baseName} \
      --remove problematic_samples.txt \
      --make-bed \
      --out B3
    '''
}

process plot_sex {
    input:
    path(sexcheck) 

    output:
    path "sexcheck.png", emit: figure

    shell:
    '''
    plot_sex.R !{sexcheck}
    '''
}

// STEP B4: Remove sex chromosomes ---------------------------------------------

process extract_autosomal {
    input:
    path(B3_bed)   
    path(B3_bim)
    path(B3_fam)
    
    output:
    path "B4.bed", emit: bed
    path "B4.bim", emit: bim 
    path "B4.fam", emit: fam
    
    shell:
    '''
    # Extract only autosomal chromosomes (for studies that don't want to 
    # include sex chromosomes)
    awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' !{B3_bim} > \
      autosomal_SNPs.txt 
    plink --bfile !{B3_bed.baseName} \
      --extract autosomal_SNPs.txt \
      --make-bed \
      --out B4    
    ''' 
}

// STEP B5: Remove SNPs with extreme heterozygosity ----------------------------
// Extract highly independent SNPs based on LD and remove MHC region 

process heterozygosity_rate {
    input:
    path(B4_bed)
    path(B4_bim)
    path(B4_fam) 
    path(exclude_regions)
    
    output:
    path "only_indep_snps.het", emit: het
    path "independent_SNPs.prune.in", emit: ind_snps // into ind_SNPs, ind_SNPs_popstrat

    shell:
    '''    
    plink --bfile !{B4_bed.baseName} \
      --exclude !{exclude_regions} \
      --indep-pairwise !{params.indep_pairwise} \
      --out independent_SNPs \
      --range
    plink --bfile !{B4_bed.baseName} \
      --extract independent_SNPs.prune.in \
      --het \
      --out only_indep_snps 
    '''
}

process plot_heterozygosity { 
    input: 
    path het

    output:
    path "het_failed_samples.txt", emit: failed
    path "heterozygosity_rate.png", emit: figure

    shell:
    '''
    plot_heterozygosity.R !{het} # get outliers too
    '''
}

process heterozygosity_prune {
    input:
    path(B4_bed)
    path(B4_bim)
    path(B4_fam)
    path(het_failed)

    output:
    path "B5.bed", emit: bed
    path "B5.bim", emit: bim
    path "B5.fam", emit: fam

    shell:
    '''
    cut -f 1,2 !{het_failed} > het_failed_plink.txt
    plink --bfile !{B4_bim.baseName} \
      --make-bed \
      --remove het_failed_plink.txt \
      --out B5
    '''
}

// STEP B6: Remove relatives ---------------------------------------------------
process relatedness {
    input:
    path(B5_bed)
    path(B5_bim)
    path(B5_fam)
    path(ind_SNPs)
    path(imiss)

    output:
    path "B6.bed", emit: bed
    path "B6.bim", emit: bim
    path "B6.fam", emit: fam
    
    shell:
    '''
    plink --bfile !{B5_bed.baseName} \
      --extract !{ind_SNPs} \
      --genome \
      --min !{params.pihat} \
      --out pihat_0.125
      
    # Identify all pairs of relatives with pihat > 0.125 and exclude one of the
    # relatives of each pair, having the most missingness. Output those failing
    # samples to pihat_failed_samples.txt
    run_IBD_QC.pl !{imiss} pihat_0.125.genome !{params.pihat}

    plink --bfile !{B5_bed.baseName} \
      --remove pihat_failed_samples.txt \
      --make-bed \
      --out B6
    '''
}

// STEP B7: Remove samples with missing phenotypes -----------------------------
process missing_phenotype {
    input:
    path(B6_bed)
    path(B6_bim)
    path(B6_fam)

    output:
    path "B7.bed", emit: bed
    path "B7.bim", emit: bim
    path "B7.fam", emit: fam

    shell:
    '''
    plink --bfile !{B6_bed.baseName} \
      --prune \
      --make-bed \
      --out B7
    '''
}

// STEP B8: Missingness per variant --------------------------------------------

process mpv {
    input:
    path(B7_bed)
    path(B7_bim)
    path(B7_fam)

    output:
    path "B8.bed", emit: bed
    path "B8.bim", emit: bim
    path "B8.fam", emit: fam
    path "plink.lmiss", emit: lmiss

    shell:
    '''
    plink --bfile !{B7_bed.baseName} --missing 
    plink --bfile !{B7_bed.baseName} \
      --geno !{params.variant_geno} \
      --make-bed \
      --out B8 
    '''
}

process plot_mpv {
  publishDir "${params.results}/qc/", mode: 'copy'

  input:
  path lmiss

  output:
  path "variant_missingness.png", emit: figure
    
  shell:
  '''
  plot_variant_missingness.R !{lmiss}
  '''
}

// STEP B9: Hardy_Weinberg equilibrium (HWE) -----------------------------------

process hardy {
  input:
  path(B8_bed)
  path(B8_bim)
  path(B8_fam)
  
  output:
  path "B9.bed", emit: bed
  path "B9.bim", emit: bim
  path "B9.fam", emit: fam
  path "plink_sub.hwe", emit: sub
  path "plinkzoomhwe.hwe", emit: zoom
  
  shell: 
  '''
  plink --bfile !{B8_bed.baseName} --hardy 
  # sample 1% of SNPs
  head -n1 plink.hwe > plink_sub.hwe
  perl -ne 'print if (rand() < 0.01)' <(tail -n +2 plink.hwe) >> plink_sub.hwe
  awk '{ if ($3=="TEST" || $3=="UNAFF" && $9 <0.001) print $0 }' \
	  plink.hwe > plinkzoomhwe.hwe
  plink --bfile !{B8_bed.baseName} \
    --hwe !{params.hwe} \
    --make-bed \
    --out B9 
  '''
}

process plot_hardy {
  publishDir "${params.results}/qc/", mode: 'copy'
  
  input:
  path sub
  path zoom

  output:
  path "*.png", emit: figure
  
  shell:
  '''
  hwe.R !{sub} ""
  hwe.R !{zoom} "strongly deviating SNPs only"
  '''
}

// STEP B10: Remove low minor allele frequency (MAF) ---------------------------

process maf {  
  input:
  path(B9_bed)
  path(B9_bim)
  path(B9_fam)

  output:
  path "B10.bed", emit: bed
  path "B10.bim", emit: bim
  path "B10.fam", emit: fam
  path "MAF_check.frq", emit: frq
  
  shell:
  '''
  plink --bfile !{B9_bed.baseName} \
    --freq \
    --out MAF_check
  plink --bfile !{B9_bed.baseName} \
    --maf !{params.maf} \
    --make-bed \
    --out B10
  '''
}

process plot_maf {
  publishDir "${params.results}/qc/", mode: 'copy'
  
  input:
  path maf_frq

  output:
  path "maf.png", emit: figure
  
  shell:
  '''
  plot_maf.R !{maf_frq}
  '''
}

// STEP B11: Test missingness in case / control status -------------------------

process test_missing {
  publishDir "${params.results}/qc/bfiles", mode: 'copy'
  
  input:
  path(B10_bed)
  path(B10_bim)
  path(B10_fam)
  
  output:
  path "B11.bed", emit: bed
  path "B11.bim", emit: bim
  path "B11.fam", emit: fam
  
  shell:
  '''
  plink --bfile !{B10_bed.baseName} \
    --test-missing
  awk '{ if ($5 < !{params.missingness}) print $2 }' plink.missing > fail_missingness.txt
  plink --bfile !{B10_bed.baseName} \
    --exclude fail_missingness.txt \
    --make-bed \
    --out B11
  '''
}
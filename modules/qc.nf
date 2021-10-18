// qc processes

// STEP C1: Remove SNPs < 90% missingness --------------------------------------
process variant_missingness {
  label 'plink1'
  
  input:
  path(in_bed)
  path(in_bim)
  path(in_fam)

  output:
  path "C1.bed", emit: bed
  path "C1.bim", emit: bim
  path "C1.fam", emit: fam
  path "C1.log", emit: log

  shell:
  '''
  plink --bfile !{in_bed.baseName} \
      --fam !{in_fam} \
      --make-bed \
      --out data 
  plink --bfile data \
      --geno 0.1 \
      --make-bed  \
      --out C1 
  '''
}

// STEP C2: Check missingness rate ---------------------------------------------
process individual_missingness {
  label 'plink1'
  
  input:
  path(C1_bed)
  path(C1_bim)
  path(C1_fam)

  output:
  path "C2.bed", emit: bed
  path "C2.bim", emit: bim
  path "C2.fam", emit: fam
  path "C2.log", emit: log
  path "before.imiss", emit: imiss_before
  path "after.imiss", emit: imiss_after
  
  
  shell:
  '''
  plink --bfile !{C1_bed.baseName} \
    --missing \
    --out before
  plink --bfile !{C1_bed.baseName} \
    --make-bed \
    --mind !{params.mind} \
    --out C2 
  plink --bfile C2 \
    --missing \
    --out after
  '''
}

process plot_missingness {
  label 'small'
  publishDir "${params.results}/qc/figures", mode: 'copy'

  input:
  path(before_imiss)
  path(after_imiss)
  val(threshold)

  output:
  path "*.png", emit: figure

  shell:
  '''
  plot_sample_missingness.R !{before_imiss} !{after_imiss} !{threshold} 
  '''
}

// STEP C3: Remove samples with sex mismatch -----------------------------------
// --check-sex requires at least one X chromosome so it has to be completed
// before excluding non-automosomal SNPs
process check_sex {
    label 'plink1'
	
    input:
    path(C2_bed)
    path(C2_bim)
    path(C2_fam)
    
    output:
    path "C3.bed", emit: bed 
    path "C3.bim", emit: bim
    path "C3.fam", emit: fam
    path "C3.log", emit: log
    path "before.sexcheck", emit: sexcheck_before
    path "after.sexcheck", emit: sexcheck_after
 
    
    shell:
    '''
    plink --bfile !{C2_bed.baseName} \
        --check-sex \
        --out before
    # Identify the samples with sex discrepancies 
    grep "PROBLEM" before.sexcheck | awk '{print $1,$2}'> \
        problematic_samples.txt
    # Delete all problematic samples
    plink --bfile !{C2_bed.baseName} \
        --remove problematic_samples.txt \
        --make-bed \
        --out C3
    plink --bfile C3 \
        --check-sex \
	--out after 
    '''
}

process plot_sex {
    label 'small'
    publishDir "${params.results}/qc/figures/", mode: 'copy'

    input:
    path(sexcheck_before)
    path(sexcheck_after) 

    output:
    path "*.png", emit: figure

    shell:
    '''
    plot_sex.R !{sexcheck_before} !{sexcheck_after}
    '''
}

// STEP C4: Remove sex chromosomes ---------------------------------------------
process extract_autosomal {
    label 'plink2'
	
    input:
    path(C3_bed)   
    path(C3_bim)
    path(C3_fam)
    
    output:
    path "C4.bed", emit: bed
    path "C4.bim", emit: bim 
    path "C4.fam", emit: fam
    path "C4.log", emit: log
    
    shell:
    '''
    # Extract only autosomal chromosomes 
    plink2 --bfile !{C3_bed.baseName} \
        --autosome \
	--make-bed \
	--out C4    
    ''' 
}

// STEP C5: Remove SNPs with extreme heterozygosity ----------------------------
process heterozygosity_rate {
    label 'plink1'
	
    input:
    path(C4_bed)
    path(C4_bim)
    path(C4_fam)
    
    output:
    path "C5_het.het", emit: het
	
    shell:
    '''    
    plink --bfile !{C4_bed.baseName} \
        --het \
        --out C5_het
    '''
}

process filter_het {
    label 'small'
	
    input:
    path(het)

    output:
    path "het_failed_samples.txt", emit: failed
    
    shell:
    '''
    filter_het.R !{het}
    '''
}

process plot_heterozygosity {
    label 'small'
    publishDir "${params.results}/qc/figures/", mode: 'copy'

    input: 
    path(het_before)
    path(het_after)

    output:
    path "*.png", emit: figure

    shell:
    '''
    plot_heterozygosity.R !{het_before} !{het_after}
    '''
}

process heterozygosity_prune {
    label 'plink1'
	
    input:
    path(C4_bed)
    path(C4_bim)
    path(C4_fam)
    path(het_failed)

    output:
    path "C5.bed", emit: bed
    path "C5.bim", emit: bim
    path "C5.fam", emit: fam
    path "C5.log", emit: log
    path "C5_after.het", emit: het

    shell:
    '''
	if !{params.heterozygosity};
	then
		cut -f 1,2 !{het_failed} > het_failed_plink.txt
		plink --bfile !{C4_bim.baseName} \
		  --make-bed \
		  --remove het_failed_plink.txt \
		  --out C5
		plink --bfile C5 \
		  --het \
		  --out C5_after
	else 
		plink --bfile !{C4_bim.baseName} \
		  --make-bed \
		  --out C5
		plink --bfile C5 \
		  --het \
		  --out C5_after
	fi
    '''
}

// STEP B6: Check for cryptic relatedness -----------------------------------
process relatedness {
    label 'plink2'  

    input:
    path(C5_bed)
    path(C5_bim)
    path(C5_fam)

    output:
    path "C6.bed", emit: bed
    path "C6.bim", emit: bim
    path "C6.fam", emit: fam
    path "C6.log", emit: log
    
    shell:
    '''
    plink2 --bfile !{C5_bed.baseName} \
        --king-cutoff !{params.king_cutoff} \
        --make-bed \
        --out C6
    '''
}

// STEP C7: Remove samples with missing phenotypes -----------------------------
process missing_phenotype {
    label 'plink1'
 
    input:
    path(C6_bed)
    path(C6_bim)
    path(C6_fam)

    output:
    path "C7.bed", emit: bed
    path "C7.bim", emit: bim
    path "C7.fam", emit: fam
    path "C7.log", emit: log

    shell:
    '''
    plink --bfile !{C6_bed.baseName} \
        --prune \
        --make-bed \
        --out C7
    '''
}

// STEP E8: Check missingness per variant --------------------------------------------
process mpv {
    label 'plink1'
	
    input:
    path(E7_bed)
    path(E7_bim)
    path(E7_fam)

    output:
    path "E8.bed", emit: bed
    path "E8.bim", emit: bim
    path "E8.fam", emit: fam
    path "E8.log", emit: log
    path "E8_before.lmiss", emit: lmiss_before
    path "E8_after.lmiss", emit: lmiss_after
	
    shell:
    '''
    plink --bfile !{E7_bed.baseName} \
        --missing \
	--out E8_before
    plink --bfile !{E7_bed.baseName} \
        --geno !{params.variant_geno} \
    	--make-bed \
	--out E8 
    plink --bfile E8 \
	--missing \
	--out E8_after
    '''
}

process plot_mpv {
    label 'small'
    publishDir "${params.results}/qc/figures/", mode: 'copy'

    input:
    path(lmiss_before)
    path(lmiss_after)
    val(threshold)

    output:
    path "*.png", emit: figure
    
    shell:
    '''
    plot_variant_missingness.R !{lmiss_before} !{lmiss_after} !{threshold} 
    '''
}

// STEP E9: Check deviation from Hardy_Weinberg equilibrium (HWE) ----------------------------
process hardy {
    label 'plink1'
	
    input:
    path(E8_bed)
    path(E8_bim)
    path(E8_fam)

    output:
    path "E9.bed", emit: bed
    path "E9.bim", emit: bim
    path "E9.fam", emit: fam
    path "E9.log", emit: log
    path "plink_sub_before.hwe", emit: sub_before
    path "plinkzoomhwe_before.hwe", emit: zoom_before
    path "plink_sub_after.hwe", emit: sub_after
    path "plinkzoomhwe_after.hwe", emit: zoom_after

    shell: 
    '''
    plink --bfile !{E8_bed.baseName} \
	--hardy \
	--out plink_before
    # sample 1% of SNPs
    head -n1 plink_before.hwe > plink_sub_before.hwe
    awk 'BEGIN {srand()} !/^$/ { if (rand() <= .01) print $0}' < plink_before.hwe >> plink_sub_before.hwe
    awk '{ if ($3=="TEST" || $3=="UNAFF" && $9 <0.001) print $0 }' \
	plink_before.hwe > plinkzoomhwe_before.hwe
    plink --bfile !{E8_bed.baseName} \
        --hwe !{params.hwe} \
        --make-bed \
        --out E9 
    plink --bfile E9 \
        --hardy \
        --out plink_after
    # sample 1% of SNPs
    head -n1 plink_after.hwe > plink_sub_after.hwe
    awk 'BEGIN {srand()} !/^$/ { if (rand() <= .01) print $0}' < plink_after.hwe >> plink_sub_after.hwe
    awk '{ if ($3=="TEST" || $3=="UNAFF" && $9 <0.001) print $0 }' \
        plink_after.hwe > plinkzoomhwe_after.hwe
    '''
}

process plot_hardy {
    label 'small'
    publishDir "${params.results}/qc/figures/", mode: 'copy'

    input:
    path(sub_before)
    path(zoom_before)
    path(sub_after)
    path(zoom_after)
    val(threshold)

    output:
    path "*.png", optional: true, emit: figure

    shell:
    '''
    plot_hwe.R !{sub_before} !{sub_after} !{threshold} "" 
    plot_hwe.R !{zoom_before} !{zoom_after} !{threshold} "strongly deviating SNPs only"
    '''
}

// STEP E10: Remove low minor allele frequency (MAF) ---------------------------
process maf {  
    label 'plink1'
    if (params.linear == true){
       publishDir "${params.results}/qc/bfiles/", pattern: "E10.*",  mode: 'copy'
    }

    input:
    path(E9_bed)
    path(E9_bim)
    path(E9_fam)

    output:
    path "E10.bed", emit: bed
    path "E10.bim", emit: bim
    path "E10.fam", emit: fam
    path "MAF_check_before.frq", emit: before
    path "MAF_check_after.frq", emit: after
    path "E10.log", emit: log

    shell:
    '''
    plink --bfile !{E9_bed.baseName} \
        --freq \
        --out MAF_check_before
    plink --bfile !{E9_bed.baseName} \
        --maf !{params.maf} \
        --make-bed \
        --out E10
    plink --bfile E10 \
        --freq \
        --out MAF_check_after
    '''
}

process plot_maf {
    label 'small'
    publishDir "${params.results}/qc/figures/", mode: 'copy'
   
    input:
    path(maf_before)
    path(maf_after)
    val(threshold)

    output:
    path "*.png", emit: figure
  
    shell:
    '''
    plot_maf.R !{maf_before} !{maf_after} !{threshold} 
    '''
}

// STEP E11: Check missingness in case / control status -------------------------
process test_missing {
    label 'plink1'
    publishDir "${params.results}/qc/bfiles/", pattern: "E11.*",  mode: 'copy'

    input:
    path(E10_bed)
    path(E10_bim)
    path(E10_fam)

    output:
    path "E11.bed", emit: bed
    path "E11.bim", emit: bim
    path "E11.fam", emit: fam
    path "before.missing", emit: before
    path "after.missing", emit: after
    path "E11.log", emit: log

    shell:
    '''
    plink --bfile !{E10_bed.baseName} \
        --test-missing \
        --out before
    awk '{ if ($5 < !{params.missingness}) print $2 }' before.missing > fail_missingness.txt
    plink --bfile !{E10_bed.baseName} \
        --exclude fail_missingness.txt \
	--make-bed \
	--out E11
    plink --bfile E11 \
        --test-missing \
        --out after
    '''
}

process plot_missing_by_cohort {
    label 'small'
    publishDir "${params.results}/qc/figures/", mode: 'copy'

    input:
    path(missing_before)
    path(missing_after)
    val(threshold)

    output:
    path "*.png", optional: true, emit: figure

    shell:
    '''
    plot_missing_cohort.R !{missing_before} !{missing_after} !{threshold} 
    ''' 
}

// STEP E12: Create covariates and plot PCA -------------------------
process pca {
    label 'plink1'
	
    input:
    path(bed)
    path(bim)
    path(fam)
    path(exclude_regions)

    output:
    path "E12_indep.log", emit: log 
    path "E12_pca.eigenvec", emit: eigenvec_user 
    path "E12_indep.fam", emit: fam
    
    shell:
    '''
    plink --bfile !{bed.baseName} \
	--exclude !{exclude_regions} \
	--indep-pairwise !{params.indep_pairwise} \
	--out indepSNPs_1k_1
    plink --bfile !{bed.baseName} \
	--extract indepSNPs_1k_1.prune.in \
	--make-bed \
	--out E12_indep
    # Perform a PCA on user's data   
    plink --bfile E12_indep \
        --pca header \
        --out E12_pca
    '''
}

process pca_covariates {
    label 'small'
    
    input:
    path(eigenvec_user)

    output:
    path "covar_pca", emit: covar
    
    shell:
    '''
    # Create covariate file including the first X PCs that the user requested
    awk -v var="!{params.pca_covars}" '{for(i=1;i<=var+2;i++) printf $i" "; print ""}' E12_pca.eigenvec > covar_pca
    '''
}

process plot_pca_user_data {
    label 'small'
    publishDir "${params.results}/qc/figures", mode: 'copy'
    
    input:
    path(eigenvec_user)
    path(fam)

    output:
    path "*.png", emit: figure
    path "*.rds", emit: rds
    
    shell:
    '''
    # Create case/control file
    awk '{print $1, $2, $6}' !{fam} > status
    plot_pca_OnlyUsersData.r !{eigenvec_user} status
    '''    
}

// STEP E13: Parse all log files and combine them into a .txt file -----------------
process parse_logs {
    label 'small'

    publishDir "${params.results}/${dir}/figures", mode: 'copy', pattern: "*.png"
    publishDir "${params.results}/${dir}/logs", mode: 'copy', pattern: "*.txt"

    input:
    val(dir)
    path(logs)
    val(fn)

    output:
    path "${fn}", emit: log
    path "*.png", emit: figure

    shell:
    '''
    echo "stage variants samples pheno pheno_case pheno_control pheno_miss wd" > !{fn}
    ls *.log | sort -V | xargs -n1 parse_logs.awk >> !{fn}
    plot_logs.R !{fn} $(basename -s .txt !{fn})
    '''
}

// STEP E14: Make an .html report ---------------------------------------------------
process report {
    label 'small'

    publishDir "${params.results}/${dir}/", mode: 'copy'
    // rmarkdown doesn't respect symlinks
    // https://github.com/rstudio/rmarkdown/issues/1508
    stageInMode 'copy'

    input:
    val(dir)
    path(x)
    path(rmd)

    output:
    path "*_report.html"

    shell:
    '''
    #!/usr/bin/env Rscript

    rmarkdown::render('!{rmd}', output_options=list(self_contained=TRUE))
    '''
}

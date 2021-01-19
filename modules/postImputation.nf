// Post-imputation (IN PROGRESS)
// ============================================================================

params.inplink = "$baseDir/../../results/sample_qc/sample_variant*"
params.inimp = "$baseDir/../../results/imputation/*"
params.outdir = "$baseDir/../../results"
outdir = params.outdir + '/postImputation/'

log.info """\
         snpQT step04: Post-imputation QC 
         """
         .stripIndent()

// User's imputed chromosomes --------------------------------------------------
// extract chromosome number from file name to make a tuple of [chr, file path]

Channel
  .fromPath( params.inimp )
  .collect()
  .set { inimp }

Channel
  .fromPath( params.inplink )
  .collect()
  .set { inplink }

// STEP E1: Merge all imputed chromosomes with bcftools, so that multi-allelics can be merged, -n is used since files are already sorted after imputation
process merge_imputed_chrom {
    input:
    file '*' from inimp

    output:
    file 'merged_imputed.vcf.gz' into E1
    
    shell:
    '''
    # file order is important so use command substition
    bcftools concat -n $(ls *.vcf.gz | sort -V) -Oz -o merged_imputed.vcf.gz
    '''
}

// STEP E2: Filter all poorly imputed variants based on info score (check impute5 output), filter based on MAF, annotate missing SNP ids, 

process filter_imputed_chrom {
    input:
    file E1

    output:
    file 'E2.bed' into E2_bed
    file 'E2.bim' into E2_bim
    file 'E2.fam' into E2_fam
    
    shell:
    '''
    plink2 --vcf !{E1} \
        --extract-if-info INFO '>'= 0.7 \
	--id-delim _ \
	--maf 0.01 \
	--set-missing-var-ids @:#:\\$r:\\$a \
	--new-id-max-allele-len 100 \
	--make-bed \
	--out E2
    '''
}

// STEP E3: Handle all categories of duplicates

process duplicates_cat1 {
    input:
    file E2_bed
    file E2_bim
    file E2_fam

    output:
    file 'E3_cat1.bed' into E3_cat1_bed
    file 'E3_cat1.bim' into E3_cat1_bim
    file 'E3_cat1.fam' into E3_cat1_fam
    
    shell:
    '''
    # Annotate all variants to this format chr:pos:ref:alt and remove exact duplicates
    plink2 --bfile !{E2_bed.baseName} \
        --set-all-var-ids @:#:\\$r:\\$a \
	--new-id-max-allele-len 1000 \
	--rm-dup force-first list \
	--make-bed \
	--out E3
    # Recover the rs ids 
    plink2 --bfile E3 \
        --recover-var-ids !{E2_bim} \
	--make-bed \
	--out E3_cat1
    '''
}

process duplicates_cat2 {
    input:
    file E3_cat1_bed
    file E3_cat1_bim
    file E3_cat1_fam

    output:
    file 'E3_cat2.bed' into E3_cat2_bed
    file 'E3_cat2.bim' into E3_cat2_bim
    file 'E3_cat2.fam' into E3_cat2_fam

    shell:
    '''
    # Identify the multi-allelics based on position and reference allele
    cut -f 1,4,6 !{E3_cat1_bim} | sort | uniq -d | cut -f 2 | grep -w -F -f - !{E3_cat1_bim} | cut -f 2 > multi_allelics.txt
    plink2 --bfile !{E3_cat1_bed.baseName} \
        --exclude multi_allelics.txt \
	--make-bed \
	--out E3_cat2
    '''
}

process duplicates_cat3 {
    input:
    file E3_cat2_bed
    file E3_cat2_bim
    file E3_cat2_fam

    output:
    file 'E3_cat3.bed' into E3_cat3_bed
    file 'E3_cat3.bim' into E3_cat3_bim
    file 'E3_cat3.fam' into E3_cat3_fam
    
    shell:
    '''
    cut -f 2 !{E3_cat2_bim} | sort | uniq -d > merged_variants.txt
    plink2 --bfile !{E3_cat2_bim.baseName} \
        -extract merged_variants.txt \
	--make-bed \
	--out merged_snps
    plink2 --bfile !{E3_cat2_bim.baseName} \
        --exclude merged_variants.txt \
	--make-bed \
	--out excluded_snps
    plink2 --bfile merged_snps \
        --set-all-var-ids @:#:\\$r:\\$a \
	--new-id-max-allele-len 1000 \
	--make-bed \
	--out annotated
    plink --bfile excluded_snps \
        --bmerge annotated \
	--make-bed \
	--out E3_cat3
    '''
}

// STEP E4: update phenotype information

process update_phenotype {
    publishDir outdir, mode: 'copy', overwrite: true
    
    input:
    file E3_cat3_bed
    file E3_cat3_bim
    file E3_cat3_fam
    file inplink

    output:
    file 'post_annotation*'
    
    shell:
    '''
    plink2 --bfile !{E3_cat3_bed.baseName} \
        --fam !{inplink[0]} \
	--make-bed \
	--out post_annotation
    '''
}
// STEP E1: Merge all imputed chromosomes with bcftools, so that multi-allelics can be merged, -n is used since files are already sorted after imputation
process merge_imp {
    input:
    path(imp)
    
    output:
    path 'merged_imputed.vcf.gz', emit: vcf
    
    shell:
    '''
    # file order is important so use command substition
    bcftools concat -n $(ls *.vcf.gz | sort -V) -Oz -o merged_imputed.vcf.gz
    '''
}

// STEP E2: Filter all poorly imputed variants based on info score (check impute5 output), filter based on MAF, annotate missing SNP ids, 

process filter_imp {
    input:
    path(imp)

    output:
    path "E2.bed", emit: bed
    path "E2.bim", emit: bim
    path "E2.fam", emit: fam
    
    shell:
    '''
    plink2 --vcf !{imp} \
        --extract-if-info INFO '>'= !{params.info} \
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
    path(bed)
    path(bim)
    path(fam)

    output:
    path "E3_cat1.bed", emit: bed
    path "E3_cat1.bim", emit: bim
    path "E3_cat1.fam", emit: fam
    
    shell:
    '''
    # Annotate all variants to this format chr:pos:ref:alt and remove exact duplicates
    plink2 --bfile !{bed.baseName} \
        --set-all-var-ids @:#:\\$r:\\$a \
	--new-id-max-allele-len 1000 \
	--rm-dup force-first list \
	--make-bed \
	--out E3
    # Recover the rs ids 
    plink2 --bfile E3 \
        --recover-var-ids !{bim} \
	--make-bed \
	--out E3_cat1
    '''
}

process duplicates_cat2 {
    input:
    path(bed)
    path(bim)
    path(fam)
    
    output:
    path "E3_cat2.bed", emit: bed
    path "E3_cat2.bim", emit: bim 
    path "E3_cat2.fam", emit: fam

    shell:
    '''
    # Identify the multi-allelics based on position and reference allele
    cut -f 1,4,6 !{bim} | sort | uniq -d | cut -f 2 | grep -w -F -f - !{bim} | cut -f 2 > multi_allelics.txt
    plink2 --bfile !{bed.baseName} \
        --exclude multi_allelics.txt \
	--make-bed \
	--out E3_cat2
    '''
}

process duplicates_cat3 {
    input:
    path(bed)
    path(bim)
    path(fam)

    output:
    path "E3_cat3.bed", emit: bed
    path "E3_cat3.bim", emit: bim 
    path "E3_cat3.fam", emit: fam
    
    shell:
    '''
    cut -f 2 !{bim} | sort | uniq -d > merged_variants.txt
    plink2 --bfile !{bim.baseName} \
        -extract merged_variants.txt \
	--make-bed \
	--out merged_snps
    plink2 --bfile !{bim.baseName} \
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
    input:
    path(bed)
    path(bim)
    path(fam)
    path(user_fam)

    output:
    path "post_annotation.bed", emit: bed
    path "post_annotation.bim", emit: bim
    path "post_annotation.fam", emit: fam
    
    shell:
    '''
    plink2 --bfile !{bed.baseName} \
        --fam !{user_fam} \
	--make-bed \
	--out post_annotation
    '''
}
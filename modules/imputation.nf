// Pre-imputation 
// =============================================================================

// STEP F1: Set chromosome codes ---------------------------------------------
process set_chrom_code {
    label 'plink2'
	
    input:
    path(bed)
    path(bim)
    path(fam)

    output:
    path "F1.bed", emit: bed
    path "F1.bim", emit: bim
    path "F1.fam", emit: fam
    path "F1.log", emit: log
    
    shell:
    '''
    plink2 --bfile !{bed.baseName} \
        --output-chr MT \
	--make-bed \
        --out F1  
    '''
}

// STEP F2: D2: Remove ambiguous SNPs and flip reverse SNPs -------------------------
// note: taken care of by popstrat modules D4 now

// STEP F3: Remove one of each pair of duplicated SNPs 
process fix_duplicates {
    label 'plink2'
	
    input:
    path(bed)
    path(bim)
    path(fam)

    output:
    path "F3.bed", emit: bed
    path "F3.bim", emit: bim
    path "F3.fam", emit: fam
    path "F3.log", emit: log
    
    shell:
    '''
    # F3: Deduplicate variants
    plink2 --bfile !{bed.baseName} \
      --rm-dup 'force-first' \
      --make-bed \
      --out F3
    '''
}

// STEP F4: Convert Plink file to .bcf file ---------------------------------
process to_vcf {
    label 'plink2'
	
    input:
    path(bed)
    path(bim)
    path(fam)
    
    output:
    path "F4.vcf.gz", emit: vcf
    path "F4.log", emit: log

    shell:
    '''
    plink2 --bfile !{bed.baseName} \
        --export vcf bgz \
        --out F4
    '''
}

// STEP F5: Convert Plink file to .bcf file ---------------------------------
process to_bcf {
    label 'bcftools'
	
    input:
    path(vcf)
    
    output:
    path "F5.bcf", emit: bcf

    shell:
    '''
    bcftools index !{vcf}
    bcftools convert !{vcf} -Ou -o F5.bcf --threads !{task.cpus}
    '''
}

// STEP F6: Check and fix the REF allele --------------------------------------
process check_ref_allele {
    label 'bcftools'
	
    input:
    path(bcf)
    path(dbsnp)
    path(dbsnp_idx)
    path(g37)

    output:
    path "F6.bcf", emit: bcf

    shell:
    '''
    bcftools +fixref !{bcf} \
        -Ob -o F6.bcf -- \
        -d -f !{g37} \
        -i !{dbsnp}
    '''
}

// STEP F7: Sort BCF, convert .bcf file to .vcf.gz file and index the vcf.gz -------------------------------
process bcf_to_vcf {
    label 'bcftools'
    publishDir "${params.results}/preImputation/files", mode: 'copy'

    input:
    path(bcf)
    
    output:
    path "F7.vcf.gz", emit: vcf
    path "F7.vcf.gz.csi", emit: idx
	
    shell:
    '''
    bcftools sort !{bcf} | bcftools convert -Oz --threads !{task.cpus} > F7.vcf.gz
    bcftools index F7.vcf.gz 
    '''
}

// Imputation 
// =============================================================================

// STEP G1: Split vcf.gz file in chromosomes and index all chroms ---------------------------------
process split_user_chrom {
    label 'bcftools'
	
    input:
    path(vcf)
    path(idx)
    each chr
    
    output:
    tuple val(chr), file('G1.vcf.gz'), file('G1.vcf.gz.csi'), emit: chrom 

    shell:
    '''
    bcftools view -r !{chr} !{vcf} -Oz -o G1.vcf.gz --threads !{task.cpus}
    bcftools index G1.vcf.gz --threads !{task.cpus}
    '''
}

// STEP G2: Perform phasing using shapeit4 --------------------------
process phasing {
    label 'shapeit'
	
    input:
    tuple val(chr), file('G1.vcf.gz'), file('G1.vcf.gz.csi'), \
        file('genetic_maps.b37.tar.gz')  

    output:
    tuple val(chr), file('G2.vcf.gz'), emit: chrom, optional: true

    shell:
    '''
    gunzip genetic_maps.b37.tar.gz
    tar -xf genetic_maps.b37.tar
    gunzip chr!{chr}.b37.gmap.gz # decompress the chromosome we need 

    # || true allows optional output without an error
    # people might not always be imputing every chromosome
    shapeit4 --input G1.vcf.gz \
        --map chr!{chr}.b37.gmap \
        --region !{chr} \
        --thread !{task.cpus} \
        --output G2.vcf.gz \
        --log log_chr.txt || true     
    '''
}

// STEP G3: Index phased chromosomes --------------------------
process bcftools_index_chr {
    label 'bcftools'
	
    input:
    tuple val(chr), path('chr.vcf.gz')

    output:
    tuple val(chr), path('chr.vcf.gz'), path('chr.vcf.gz.csi'), emit: chrom_idx

    shell:
    '''
    bcftools index chr.vcf.gz 
    '''
}

// STEP G4: Tabix reference files ----------------------------
process tabix_chr {
    label 'small'
    label 'tabix'
	
	
    input:
    tuple val(chr), path('chr.vcf.gz')

    output:
    tuple val(chr), path('chr.vcf.gz'), path('chr.vcf.gz.tbi'), emit: chrom_idx

    shell:
    '''
    tabix -p vcf chr.vcf.gz
    '''
}

// STEP G5: Convert vcf reference genome into a .imp5 format for each chromosome
process convert_imp5 {
    label 'impute5'
	
    input:
    tuple val(chr), file('ref_chr.vcf.gz'), file('ref_chr.vcf.gz.tbi')

    output:
    tuple val(chr), file('1k_b37_reference_chr.imp5'), \
        file('1k_b37_reference_chr.imp5.idx'), emit: chrom

    shell:
    '''
    imp5Converter!{params.impute5_version} --h ref_chr.vcf.gz \
        --r !{chr} \
	--threads !{task.cpus} \
        --o 1k_b37_reference_chr.imp5
    '''
}

// STEP G6: Perform imputation using impute5 ---------------------------------
// join phased vcfs with imp5 based on chrom value 
// then combine so each tuple element has a shapeit4 map file 

process impute5 {
    label 'bigmem'
    label 'impute5'
	
    input:
    tuple val(chr), file('1k_b37_reference_chr.imp5'), \
        file('1k_b37_reference_chr.imp5.idx'), file('G2.vcf.gz'), \
        file('G2.vcf.gz.csi'), file('genetic_maps.b37.tar.gz') 
       
    output:
    path "${chr}.vcf.gz", emit: imputed

    shell:
    '''
    tar -xzf genetic_maps.b37.tar.gz
    gunzip chr!{chr}.b37.gmap.gz # decompress the chromosome we need

    impute5!{params.impute5_version} --h 1k_b37_reference_chr.imp5 \
        --m chr!{chr}.b37.gmap \
        --g G2.vcf.gz \
        --r !{chr} \
        --out-gp-field \
        --o !{chr}.vcf.gz
    '''
}

// STEP G7: Merge all imputed chromosomes with bcftools, so that multi-allelics can be merged, -n is used since files are already sorted after imputation
process merge_imp {
    label 'bcftools'
	
    input:
    path(imp)
    
    output:
    path 'merged_imputed.vcf.gz', emit: vcf
    
    shell:
    '''
    # file order is important so use command substition
    bcftools concat -n $(ls *.vcf.gz | sort -t . -k 1n) -Oz -o merged_imputed.vcf.gz --threads !{task.cpus}
    '''
}

// Finished!


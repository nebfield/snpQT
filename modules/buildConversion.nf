// Step B1: Decompressing fasta file (workflow download_db.nf)-----------------
// Step B2: Create a dictionary file ------------------------------------------
process dictionary {
    label 'small'
    
    input:
    path(fa)

    output:
    path "*.dict", emit: dict

    shell:
    '''
    # make it so!
    picard \
	CreateSequenceDictionary \
	-R !{fa} \
	-O !{fa.baseName}.fa.dict
    '''
}

// STEP B3: Change the chr ids ------------------------------------------------
process num_to_chr {
    input:
    path(in_vcf)
    path(chr_map)

    output:
    path "out.vcf.gz", emit: vcf

    shell:
    '''
    bcftools annotate --rename-chr !{chr_map} !{in_vcf} \
	-Oz -o out.vcf.gz --threads !{task.cpus}
    '''
}

// STEP B4: Run liftOver to map genome build -----------------------------------
process liftover {
    label 'small'
    memory { (params.mem * 1.2) + 'G' }

    input:
    path(vcf)
    path(hg)
    path(chain)
    path(dict)

    output:
    path "out.vcf", emit: vcf
    path "rejected_variants.vcf", emit: rejected_vcf

    shell:
    '''
    # !{dict} unused but needed to stage in file
    picard "-Xmx!{params.mem}G" LiftoverVcf \
	-I !{vcf} \
	-O out.vcf \
	-CHAIN !{chain} \
	-REJECT rejected_variants.vcf \
	-R !{hg}
    '''
}

// STEP B5: Reverse Chr1To1 ---------------------------------------------------
process chr_to_num {  
    publishDir "${params.results}/convertBuild/files/", mode: 'copy'

    input:
    path(vcf)
    path(chr_map)

    output:
    path "converted.vcf.gz", emit: vcf

    shell:
    '''
    awk '{print $2 "\t" $1}' !{chr_map} > Chr1To1.txt

    # Change the chromosome ids again
    bcftools annotate --rename-chr Chr1To1.txt !{vcf} -Oz -o converted.vcf.gz --threads !{task.cpus}
    '''
}

// STEP B6: Convert VCF to PLINK format ---------------------------------------
process vcf_to_plink {
    publishDir "${params.results}/convertBuild/files/", mode: 'copy'

    input:
    path(vcf)
    path(in_fam)

    output:
    path "converted.bed", emit: bed
    path "converted.bim", emit: bim
    path "converted.fam", emit: fam

    shell:
    '''
    # Convert VCF to PLINK format
    plink2 --vcf !{vcf} \
        --id-delim _ \
        --max-alleles 2 \
        --chr 1-22 XY \
        --allow-extra-chr \
        --make-bed \
        --out converted

    # Update phenotypes
    plink2 --bfile converted \
        --fam !{in_fam} \
        --make-bed \
        --out converted 
    '''
}

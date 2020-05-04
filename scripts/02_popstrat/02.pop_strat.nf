params.inbed = "../../results/sample_qc/sample_variant_qc.*"
params.inbim = "../../results/sample_qc/sample_variant_qc.bim"
params.infam = "../../results/sample_qc/sample_variant_qc.fam"
params.refbed = "$SNPQT_DB_DIR/1kG_PCA5.bed"
params.refbim = "$SNPQT_DB_DIR/1kG_PCA5.bim"
params.reffam = "$SNPQT_DB_DIR/1kG_PCA5.fam"
params.outdir = "$baseDir/../../results"
outdir = params.outdir + '/pop_strat/'

log.info """\
         snpQT step02: population stratification 
         input file: ${params.inbim}
         outdir: ${params.outdir}
         """
         .stripIndent()

Channel.fromPath( params.inbed ).set { inbed } 
Channel.fromPath( params.inbim ).set { inbim } 
Channel.fromPath( params.infam ).set { infam } 
Channel.fromPath( params.refbed ).set { refbed } 
Channel.fromPath( params.refbim ).set { refbim } 
Channel.fromPath( params.reffam ).set { reffam } 

// STEP C3: filter minor allele frequency from user's dataset ------------------

process filter_maf {
    input:
    file inbed 
    file inbim
    file infam

    output:
    file "maf_filtered*" into maf_filtered

    """
    plink --bfile sample_variant_qc --maf 0.05 --make-bed --out maf_filtered
    """
}

// STEP C4: harmonise 1000 genomes data  ---------------------------------------

process intersect_variants {
    input:
    file maf_filtered
    file refbed
    file refbim
    file reffam

    output:
    file "user_intersected*" into user_intersected
    file "1kG_PCA5_intersected*" into ref_intersected

    shell:
    '''
    awk '{print $2}' maf_filtered.bim > user_snps.txt
    plink --bfile 1kG_PCA5 --extract user_snps.txt --make-bed \
      --out 1kG_PCA5_intersected
    # Extract the variants present in 1000 Genomes dataset from the  dataset.
    awk '{print $2}' 1kG_PCA5_intersected.bim > 1kG_PCA5_SNPs.txt
    plink --bfile maf_filtered --extract 1kG_PCA5_SNPs.txt --recode \
      --make-bed --out user_intersected
    # The datasets now contain the exact same variants.
    '''
}

process harmonise_build {
    input:
    file user_intersected
    file ref_intersected 
    
    shell:
    '''
    ## The datasets must have the same build. Change the build 1000 Genomes data build.
    awk '{print$2,$4}' user_intersected.map > build_map.txt
    # build_map.txt contains one SNP-id and physical position per line.
    plink --bfile 1kG_PCA5_plink_12 --update-map build_map.txt --make-bed --out 1kG_PCA5_plink_12
    # 1kG_PCA5_plink_12 and plink_12_1 now have the same build.

    # set 1k genome as reference to user's data 
    awk '{print$2,$5}' 1kG_PCA5_plink_12.bim > 1kg_ref-list.txt
    plink --bfile plink_PCA --reference-allele 1kg_ref-list.txt --make-bed --out plink_PCA-adj
    # The 1kG_PCA5_plink_12 and the plink_PCA-adj have the same reference genome for all SNPs.
    # This command will generate some warnings for impossible A1 allele assignment.
    '''
}

/* process flip_snps {
    shell:
    '''
    
#Run snpflip to identify ambiguous SNPs and SNPs that are located on the reverse strand first on user's dataset
snpflip -b plink_PCA-adj.bim -f /home/.../refs/human_g1k_v37.fasta -o plink_PCA-adj_snpflip
#outputs two files, one plink_PCA-adj_snpflip.ambiguous and plink_PCA-adj_snpflip.reverse

    # Flip all reversed SNPs
    plink --bfile  plink_PCA-adj --flip plink_PCA-adj_snpflip.reverse --reference-allele 1kg_ref-list.txt --make-bed --out plink_PCA-adj_flipped

    #Remove ambiguous SNPs
    plink --bfile  plink_PCA-adj_flipped --exclude plink_PCA-adj_snpflip.ambiguous  --make-bed --out plink_PCA_Corr

    #Run snpflip on 1k dataset
    snpflip -b 1kG_PCA5_plink_12.bim -f /home/.../refs/human_g1k_v37.fasta -o 1kG_PCA5_plink_12_snpflip

    # Flip all reversed SNPs
    plink --bfile  1kG_PCA5_plink_12 --flip 1kG_PCA5_plink_12_snpflip.reverse --reference-allele 1kg_ref-list.txt --make-bed --out 1kG_PCA5_plink_12_flipped

    #Remove ambiguous SNPs
    plink --bfile  1kG_PCA5_plink_12 --exclude plink_PCA-adj_snpflip.ambiguous  --make-bed --out 1kG_PCA5_Corr
    '''
}
 */

// STEP C5: merge 1000 genomes data  -------------------------------------------

// STEP C6: PCA anchored on 1000 genomes  --------------------------------------
// TODO

// STEP C7: make racefile  -----------------------------------------------------
// TODO

// STEP C8: plot PCA  ----------------------------------------------------------
// TODO

// STEP C9: extract homogenous ethnic group  -----------------------------------
// TODO

// STEP C10: Logistic regression  ----------------------------------------------
// TODO
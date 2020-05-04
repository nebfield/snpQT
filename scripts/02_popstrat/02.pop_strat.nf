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

Channel.fromPath("$SNPQT_DB_DIR/human_g1k_v37.fasta").set{ g37 }

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
    
    output:
    file "user_intersected-adj*" into user_harmonised
    file "1kG_PCA5_intersected_mapped*" into ref_harmonised
    file "1kg_ref-list.txt" into ref_list

    shell:
    '''
    ## The datasets must have the same build. Change the build 1000 Genomes data build.
    awk '{print$2,$4}' user_intersected.map > build_map.txt
    # build_map.txt contains one SNP-id and physical position per line.
    plink --bfile 1kG_PCA5_intersected --update-map build_map.txt \
      --make-bed --out 1kG_PCA5_intersected_mapped

    # set 1k genome as reference to user's data 
    awk '{print$2,$5}' 1kG_PCA5_intersected_mapped.bim > 1kg_ref-list.txt
    plink --bfile user_intersected --reference-allele 1kg_ref-list.txt \
      --make-bed --out user_intersected-adj
    # The 1kG_PCA5_intersected and the user_intersected-adj have the same 
    # reference genome for all SNPs
    # This command will generate some warnings for impossible A1 allele assignment.
    '''
}

process flip_snps {
    input:
    file user_harmonised
    file ref_harmonised
    file ref_list
    file g37

    output:
    file "plink_PCA_Corr_2*" into user_flipped
    file "1kG_PCA5_Corr_2*" into ref_flipped
    
    shell:
    '''
    # Run snpflip to identify ambiguous SNPs and SNPs that are located on 
    # the reverse strand first on user's dataset
    snpflip -b user_intersected-adj.bim \
      -f !{g37} \
      -o user_intersected-adj_snpflip
    # outputs two files, one user_intersected-adj_snpflip.ambiguous 
    # user_intersected-adj_snpflip.reverse

    # Flip all reversed SNPs
    plink --bfile user_intersected-adj \
      --flip user_intersected-adj_snpflip.reverse \
      --reference-allele 1kg_ref-list.txt \
      --make-bed --out plink_PCA-adj_flipped

    #Remove ambiguous SNPs
    plink --bfile plink_PCA-adj_flipped \
      --exclude user_intersected-adj_snpflip.ambiguous \
      --make-bed --out plink_PCA_Corr

    # Run snpflip on 1k dataset
    snpflip -b 1kG_PCA5_intersected_mapped.bim \
      -f !{g37} \
      -o 1kG_PCA5_intersected_snpflip

    # Flip all reversed SNPs
    plink --bfile 1kG_PCA5_intersected_mapped \
      --flip 1kG_PCA5_intersected_snpflip.reverse \
      --reference-allele 1kg_ref-list.txt \
      --make-bed --out 1kG_PCA5_intersected_flipped

    # Remove ambiguous SNPs
    plink --bfile 1kG_PCA5_intersected_flipped \
      --exclude 1kG_PCA5_intersected_snpflip.ambiguous \
      --make-bed --out 1kG_PCA5_Corr

    # Find differences between the two files that still appeat after flipping an removing ambiguous SNPs
    awk '{print$2,$5,$6}' plink_PCA_Corr.bim > user_data_corrected_tmp
    awk '{print$2,$5,$6}' 1kG_PCA5_Corr.bim > 1k_corrected_tmp
    sort user_data_corrected_tmp 1k_corrected_tmp | uniq \
      -u  > uncorresponding_SNPs.txt

    # Keep only the unique SNP ids 
    awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exclusion.txt

    # Remove the problematic SNPs from both datasets.
    plink --bfile plink_PCA_Corr \
      --exclude SNPs_for_exclusion.txt \
      --make-bed --out plink_PCA_Corr_2
    plink --bfile 1kG_PCA5_Corr \
      --exclude SNPs_for_exclusion.txt \
      --make-bed --out 1kG_PCA5_Corr_2
    '''
}
 
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
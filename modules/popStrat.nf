// STEP C3: filter minor allele frequency from user's dataset ------------------

process filter_maf {
    input:
    path(bed)
    path(bim)
    path(fam)
    path(db)
    
    output:
    path "C3.bed", emit: bed
    path "C3.bim", emit: bim
    path "C3.fam", emit: fam
    
    shell:
    '''
    # get specific files from db
    exclude=!{db}"/PCA.exclude.regions.b37.txt"
    # MAF filtering < 5%
    plink --bfile !{bed.baseName} \
      --maf 0.05 \
      --make-bed \
      --out maf_filtered
    
    # Pruning in user's dataset
    plink --bfile maf_filtered \
      --exclude $exclude \
      --indep-pairwise 50 5 0.2 \
      --out indepSNPs_1k
      
    plink --bfile maf_filtered \
      --extract indepSNPs_1k.prune.in \
      --make-bed \
      --out C3
   '''
}

// STEP C4: Fix strand errors and remove ambiguous SNPs ------------------------

process run_snpflip {
  input: 
  path(bed)
  path(bim)
  path(fam)
  path(db)
  
  output:
  path "flipped_snps.reverse", emit: rev
  path "flipped_snps.ambiguous", emit: ambig
  
  shell:
  '''
  # extract files from db
  g37=!{db}"/human_g1k_v37.fasta"
  # Run snpflip to identify ambiguous SNPs and SNPs that are located on the 
  # reverse strand first on user's dataset
  snpflip -b !{bim} \
    -f $g37 \
    -o flipped_snps
  '''
}

process flip_snps {  
  input:
  path(bed)
  path(bim)
  path(fam)
  path(snpflip_rev)
  path(snpflip_ambig)
  
  output:
  path "C4.bed", emit: bed
  path "C4.bim", emit: bim
  path "C4.fam", emit: fam
  
  shell:
  '''
  # Flip all reversed SNPs
  plink --bfile !{bed.baseName} \
    --flip !{snpflip_rev} \
    --make-bed \
    --out flipped

  # Remove ambiguous SNPs
  plink --bfile flipped \
    --exclude !{snpflip_ambig} \
    --make-bed \
    --out C4
  '''
}

// STEP C5: Align the reference allele according to 1k reference genome --------

process align {
  input:
  path(bed)
  path(bim)
  path(fam)
  path(ref_bed)
  path(ref_bim)
  path(ref_fam)

  output:
  path "C5.bed", emit: bed
  path "C5.bim", emit: bim
  path "C5.fam", emit: fam
  
  shell:
  '''
  # set 1k genome as reference to user's data 
  awk '{print$2,$5}' !{ref_bim} > 1kg_ref-list.txt

  plink --bfile !{bed.baseName} \
    --reference-allele 1kg_ref-list.txt \
    --make-bed \
    --out C5
  '''
}

// STEP C6: Merge user's dataset with 1k reference genome ----------------------

process merge {
    input:
    path(bed)
    path(bim)
    path(fam)
    path(ref_bed)
    path(ref_bim)
    path(ref_fam)
    
    output:
    path "C6.bed", emit: bed
    path "C6.bim", emit: bim
    path "C6.fam", emit: fam
    
    shell:
    '''
    # Extract the variants present in dataset from the 1000 genomes dataset
    awk '{print $2}' !{bim} > user_snps.txt
    plink --bfile !{ref_bed.baseName} \
      --extract user_snps.txt \
      --make-bed \
      --out 1kG_subset

    # Extract the variants present in 1000 Genomes dataset from the dataset
    awk '{print $2}' 1kG_subset.bim > 1kG_PCA6_SNPs.txt
    plink --bfile !{bim.baseName} \
      --extract 1kG_PCA6_SNPs.txt \
      --make-bed \
      --out C5_subset

    # The datasets now contain the exact same variants.

    # Find differences between the two files that still appeat after flipping 
    # an removing ambiguous SNPs
    awk '{print $2,$5,$6}' C5_subset.bim > user_data_corrected_tmp
    awk '{print $2,$5,$6}' 1kG_subset.bim > 1k_corrected_tmp
    sort user_data_corrected_tmp 1k_corrected_tmp | uniq -u > uncorresponding_SNPs.txt

    # Keep only the unique SNP ids 
    awk '{print $1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exclusion.txt

    # Remove the problematic SNPs from both datasets.
    plink --bfile C5_subset \
      --exclude SNPs_for_exclusion.txt \
      --make-bed \
      --out C5_subset_exclude
    plink --bfile 1kG_subset \
      --exclude SNPs_for_exclusion.txt \
      --make-bed \
      --out 1kG_subset_exclude

    # Merge user's dataset with 1000 Genomes Data
    plink --bfile C5_subset_exclude \
      --bmerge 1kG_subset_exclude.bed 1kG_subset_exclude.bim 1kG_subset_exclude.fam \
      --allow-no-sex \
      --make-bed \
      --out C6 
    '''
}

// STEP C7: PCA ----------------------------------------------------------------

process pca_prep {    
    input:
    path(bed)
    path(bim)
    path(fam)
    path(db)
    
    output:
    path("C6_indep.bed"), emit: bed 
    path("C6_indep.bim"), emit: bim
    path("C6_indep.fam"), emit: fam

    shell:
    '''
    # extract specific file from db
    exclude=!{db}"/PCA.exclude.regions.b37.txt"
    
    # recalculate independent snps
    plink --bfile !{bed.baseName} \
      --exclude $exclude \
      --indep-pairwise 50 5 0.2 \
      --out independent_SNPs 

    # Pruning on merged dataset
    plink --bfile !{bed.baseName} \
      --extract independent_SNPs.prune.in \
      --make-bed \
      --out C6_indep 
    '''
}

process racefile {
    input:
    path(db)
    
    output:
    path "super_racefile.txt", emit: super
    path "sub_racefile.txt", emit: sub

    shell:
    '''
    panel=!{db}"/integrated_call_samples_v3.20130502.ALL.panel"
    # Make 1st racefile, using the 20130502 panel using superpopulation codes
    # (i.e., AFR,AMR,EASN,SAS and EUR).
    awk '{print $1,$1,$3}' $panel > super_racefile.txt

    # Make 2nd racefile, using the 20130502 panel using subpopulation codes 
    awk '{print $1,$1,$2}' $panel > sub_racefile.txt  
    '''
}

// STEP C8: Eigensoft ----------------------------------------------------------------
process eigensoft {
    input:
    path bed
    path bim
    path fam
    path super_racefile
    path pca_fam

    output:
    path "eigenvec", emit: eigenvec
    path "merged_super_racefile.txt", emit: merged_racefile 
    path "keep_sample_list.txt", emit: keep_samples
    
    shell:
    '''
    # Concatenate racefiles: User's + super_racefile
    awk '{print $1,$2,"OWN"}' !{pca_fam} > racefile_toy_own.txt
    cat !{super_racefile} racefile_toy_own.txt | sed -e '1i\\FID IID race' >  merged_super_racefile.txt

    # Assign populations to FID and IIDs, make .pedind
    awk 'NR==FNR {h[\$2] = \$3; next} {print \$1,\$2,\$3,\$4,\$5,h[\$2]}' merged_super_racefile.txt !{fam} > C6_indep.pedind

    # make poplist.txt
    echo "OWN" > poplist.txt 
    cut -d ' ' -f 6  C6_indep.pedind | sort | uniq | grep -v 'OWN' >> poplist.txt
    
    echo "genotypename: !{bed}" > parfile
    echo "snpname:      !{bim}" >> parfile
    echo "indivname:    C6_indep.pedind" >> parfile
    echo "evecoutname:  eigenvec" >> parfile
    echo "evaloutname:  eigenval" >> parfile
    echo "numthreads:   10" >> parfile
    echo "poplistname: poplist.txt" >> parfile
    echo "numoutlierevec: 5" >> parfile
    echo "autoshrink: YES" >> parfile

    smartpca -p parfile > log.txt

    awk -F " " '{print $1}' eigenvec | sed '1d' | awk -F ":" '{print $1,$2}' > keep_sample_list.txt
    '''
}

process plot_pca {
    input:
    path eigenvec
    path super_race
    
    output:
    path "*.png"
    path "*.html"
    path "popStrat_files"

    shell:
    '''
    pop_strat.R !{eigenvec} !{super_race}
    '''
}

// STEP C9: Extract homogenous ethnic group ------------------------------------------------------

process extract_homogenous {
    input:
    path(bed)
    path(bim)
    path(fam)
    path(keep)
    
    output:
    path "C9.bed", emit: bed
    path "C9.bim", emit: bim
    path "C9.fam", emit: fam
    
    shell:
    '''
    plink -bfile !{bed.baseName} --keep !{keep} --make-bed --out C9
    '''
}

// STEP C10: Covariates ---------------------------------------------------------------------------

process pca_covariates {
    input:
    path(bed)
    path(bim)
    path(fam)
    path(db)

    output:
    path "covar_pca", emit: covar
    
    shell:
    '''
    # get specific files from db
    exclude=!{db}"/PCA.exclude.regions.b37.txt"
 
    plink --bfile !{bed.baseName} \
      --exclude $exclude \
      --indep-pairwise 50 5 0.2 \
      --out indepSNPs_1k_1
    plink --bfile !{bed.baseName} \
      --extract indepSNPs_1k_1.prune.in \
      --make-bed \
      --out C10_indep

    # Perform a PCA on user's data without ethnic outliers.    
    plink --bfile C10_indep --pca header --out C10_pca
    # Create covariate file including the first 3 PCs
    awk '{print $1, $2, $3, $4, $5}' C10_pca.eigenvec > covar_pca
    '''
}
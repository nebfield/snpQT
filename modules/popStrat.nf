// TODO: project name parameter
// TODO: hardcoded input files are called sample_variant_qc.*
// collect() and get baseName?
// TODO: big log file

params.inbed = "../../results/sample_qc/sample_variant_qc.*"
params.inbim = "../../results/sample_qc/sample_variant_qc.bim"
params.infam = "../../results/sample_qc/sample_variant_qc.fam"
params.outdir = "$baseDir/../../results"
outdir = params.outdir + '/popStrat/'

log.info """\
         snpQT step02: population stratification 
         input file: ${params.inbim}
         outdir: ${params.outdir}
         """
         .stripIndent()

Channel.fromPath( params.inbed ).set { inbed } 
Channel.fromPath( params.inbim ).set { inbim } 
Channel.fromPath( params.infam ).set { infam } 
inbed.into {inbed_maf ; inbed_extract }
inbim.into { inbim_maf ; inbim_extract }
infam.into { infam_maf ; infam_extract }
Channel.fromPath("$SNPQT_DB_DIR/all_phase3_10.bed").set { refbed }
Channel.fromPath("$SNPQT_DB_DIR/all_phase3_10.bim").set { refbim } 
Channel.fromPath("$SNPQT_DB_DIR/all_phase3_10.fam").set { reffam }
refbed.into { refbed_align ; refbed_intersect }
refbim.into { refbim_align ; refbim_intersect }
reffam.into { reffam_align ; reffam_intersect }
Channel.fromPath("$SNPQT_DB_DIR/human_g1k_v37.fasta").set{ g37 }
Channel.fromPath("$SNPQT_DB_DIR/1kG_race.txt").set{ racefile }
Channel.fromPath("$SNPQT_DB_DIR/PCA.exclude.regions.b37.txt").set{exclude} 
exclude.into { C3_exclude_regions ; C7_exclude_regions ; C10_exclude_regions }

// STEP C3: filter minor allele frequency from user's dataset ------------------

process filter_maf {
    input:
    file inbed_maf 
    file inbim_maf
    file infam_maf
    file C3_exclude_regions

    output:
    file "C3*" into C3, C3_snpflip, C3_pca
    file "C3.fam" into C3_pca_fam
    file "C3.log" into C3_log
    
    shell:
    '''
    # MAF filtering < 5%
    plink --bfile sample_variant_qc \
      --maf 0.05 \
      --make-bed \
      --out maf_filtered
    
    # Pruning in user's dataset
    plink --bfile maf_filtered \
      --exclude !{C3_exclude_regions} \
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
  container 'snpflip'

  input: 
  file C3
  file g37

  output:
  file "plink_PCA-adj_snpflip*" into snpflip_output

  shell:
  '''
  # Run snpflip to identify ambiguous SNPs and SNPs that are located on the 
  # reverse strand first on user's dataset
  snpflip -b C3.bim \
    -f !{g37} \
    -o plink_PCA-adj_snpflip
  '''
}

process flip_snps {  
  input:
  file snpflip_output
  file C3_snpflip 

  output:
  file "C4*" into C4
  file "C4.log" into C4_log
  
  shell:
  '''
  # Flip all reversed SNPs
  plink --bfile C3 \
    --flip plink_PCA-adj_snpflip.reverse \
    --make-bed \
    --out flipped

  # Remove ambiguous SNPs
  plink --bfile flipped \
    --exclude plink_PCA-adj_snpflip.ambiguous \
    --make-bed \
    --out C4
  '''
}

// STEP C5: Align the reference allele according to 1k reference genome --------

process align {
  input:
  file C4 
  file refbed_align
  file refbim_align
  file reffam_align

  output:
  file "C5*" into C5
  file "C5.log" into C5_log
  
  shell:
  '''
  # set 1k genome as reference to user's data 
  awk '{print$2,$5}' !{refbim_align} > 1kg_ref-list.txt

  plink --bfile C4 \
    --reference-allele 1kg_ref-list.txt \
    --make-bed \
    --out C5
  '''
}

// STEP C6: Merge user's dataset with 1k reference genome ----------------------

process intersect_variants {
    input:
    file C5    
    file refbed_intersect
    file refbim_intersect
    file reffam_intersect 

    output:
    file "C6*" into C6_pca, C6_racefile
    file "C6.log" into C6_log
    
    shell:
    '''
    # Extract the variants present in dataset from the 1000 genomes dataset
    awk '{print $2}' C5.bim > user_snps.txt
    plink --bfile !{refbed_intersect.baseName} \
      --extract user_snps.txt \
      --make-bed \
      --out 1kG_subset

    # Extract the variants present in 1000 Genomes dataset from the dataset
    awk '{print $2}' 1kG_subset.bim > 1kG_PCA6_SNPs.txt
    plink --bfile C5 \
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
    conda 'bioconda:eigensoft'
    
    input:
    file C6_pca
    file C7_exclude_regions 
    
    output:
    file "C6_indep.bim" into C6_indep_bim
    file "C6_indep.bed" into C6_indep_bed
    file "C6_indep.fam" into C6_indep_fam

    shell:
    '''
    # recalculate independent snps
    plink --bfile C6 \
      --exclude !{C7_exclude_regions} \
      --indep-pairwise 50 5 0.2 \
      --out independent_SNPs 

    # Pruning on merged dataset
    plink --bfile C6 \
      --extract independent_SNPs.prune.in \
      --make-bed \
      --out C6_indep 
    '''
}

process get_racefile {
    output:
    file "super_racefile.txt" into super_racefile
    file "sub_racefile.txt" into sub_racefile

    shell:
    '''
    curl -O ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
    # Make 1st racefile, using the 20130502 panel using superpopulation codes
    # (i.e., AFR,AMR,EASN,SAS and EUR).
    awk '{print $1,$1,$3}' integrated_call_samples_v3.20130502.ALL.panel > super_racefile.txt

    # Make 2nd racefile, using the 20130502 panel using subpopulation codes 
    awk '{print $1,$1,$2}' integrated_call_samples_v3.20130502.ALL.panel > sub_racefile.txt  
    '''
}

// STEP C8: Eigensoft ----------------------------------------------------------------
process eigensoft {
    container 'quay.io/biocontainers/eigensoft:7.2.1--h1d3628b_2'
 
    input:
    file C3_pca_fam
    file C6_indep_bed
    file C6_indep_bim
    file C6_indep_fam
    file super_racefile

    output:
    file "eigenvec" into C8_eigenvec
    file "merged_super_racefile.txt" into super_race_plot
    file "keep_sample_list.txt" into keep_sample_list
    
    shell:
    '''
    # Concatenate racefiles: User's + super_racefile
    awk '{print $1,$2,"OWN"}' !{C3_pca_fam} > racefile_toy_own.txt
    cat !{super_racefile} racefile_toy_own.txt | sed -e '1i\\FID IID race' >  merged_super_racefile.txt

    # Assign populations to FID and IIDs, make .pedind
    awk 'NR==FNR {h[\$2] = \$3; next} {print \$1,\$2,\$3,\$4,\$5,h[\$2]}' merged_super_racefile.txt !{C6_indep_fam} > C6_indep.pedind

    # make poplist.txt
    echo "OWN" > poplist.txt 
    cut -d ' ' -f 6  C6_indep.pedind | sort | uniq | grep -v 'OWN' >> poplist.txt
    
    echo "genotypename: !{C6_indep_bed}" > parfile
    echo "snpname:      !{C6_indep_bim}" >> parfile
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
    publishDir outdir, mode: 'copy', overwrite: true
    input:
    file C8_eigenvec
    file super_race_plot

    output:
    file "*.png"
    file "*.html"
    file "popStrat_files"

    shell:
    '''
    pop_strat.R !{C8_eigenvec} !{super_race_plot}
    '''
}

// STEP C9: Extract homogenous ethnic group ------------------------------------------------------

process extract_homogenous {
    input:
    file C3_pca
    file keep_sample_list

    output:
    file "C9*" into C9
    file "C9.log" into C9_log
    
    shell:
    '''
    plink -bfile C3 --keep !{keep_sample_list} --make-bed --out C9
    '''
}

// STEP C10: Covariates ---------------------------------------------------------------------------

process pca_covariates {
    input:
    file C9
    file C10_exclude_regions

    output:
    file "C10_pca.log" into C10_log
    
    shell:
    '''
    plink --bfile C9 --exclude !{C10_exclude_regions} --indep-pairwise 50 5 0.2 --out indepSNPs_1k_1
    plink --bfile C9 --extract indepSNPs_1k_1.prune.in --make-bed --out C10_indep
    # Perform a PCA on user's data without ethnic outliers.
    
    plink --bfile C10_indep --pca header --out C10_pca
    # Create covariate file including the first 3 PCs
    awk '{print $1, $2, $3, $4, $5}' C10_pca.eigenvec > covar_pca
    '''
}

process logs {
    input:
    file C3_log
    file C4_log
    file C5_log
    file C6_log
    file C9_log
    file C10_log

    output:
    file "popstrat.log"
    
    shell:
    '''
    ls *.log | sort -V | xargs -d '\n' grep loaded > popstrat.log
    '''
}
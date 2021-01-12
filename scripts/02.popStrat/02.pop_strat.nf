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
exclude.into { C3_exclude_regions ; C7_exclude_regions; C10_excluded_regions}

// STEP C3: filter minor allele frequency from user's dataset ------------------

process filter_maf {
    input:
    file inbed_maf 
    file inbim_maf
    file infam_maf
    file C3_exclude_regions

    output:
    file "C3*" into C3, C3_snpflip 

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

process pca {
    input:
    file C6_pca
    file C7_exclude_regions 
    
    output:
    file "C7.eigenvec" into C7_eigenvec
    file "C7*" into C7 

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
	
    # Perform PCA
    plink --bfile C6_indep \
      --pca header \
      --out C7
    '''
}

process racefile {
    input:
    file C6_racefile

    output:
    file "merged_super_racefile.txt" into super_racefile
    file "merged_sub_racefile.txt" into sub_racefile

    shell:
    '''
    curl -O ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
    # Make 1st racefile, using the 20130502 panel using superpopulation codes
    # (i.e., AFR,AMR,EASN,SAS and EUR).
    awk '{print $1,$1,$3}' integrated_call_samples_v3.20130502.ALL.panel > super_racefile_1k.txt

    # Make 2nd racefile, using the 20130502 panel using subpopulation codes 
	  awk '{print $1,$1,$2}' integrated_call_samples_v3.20130502.ALL.panel > sub_racefile_1k.txt

    # Create a racefile with user's data.
	  awk '{print$1,$2,"OWN"}' C6.fam > racefile_own.txt

	  # Concatenate racefiles: User's + super_racefile.
	  cat super_racefile_1k.txt racefile_own.txt | sed -e '1iFID IID race' > merged_super_racefile.txt

	  # Concatenate racefiles: User's + sub_racefile.
	  cat sub_racefile_1k.txt racefile_own.txt | sed -e '1iFID IID race' > merged_sub_racefile.txt    
    '''
}

// STEP C8: plot PCA  ----------------------------------------------------------

process plot_pca {
    publishDir outdir, mode: 'copy', overwrite: true, pattern: "*.png"

    input:
    file C7_eigenvec
    file super_racefile

    output:
    file "PCA.png"
    file "EUR_PCA_merge" into ethnic_cluster

    shell:
    '''
    pop_strat.R !{C7_eigenvec} !{super_racefile}

    # Exclude ethnic outliers
    # -----------------------
	  # Select individuals in plink data below cut-off thresholds. The cut-off 
    # levels are not fixed thresholds but have to be determined based on the 
    # visualization of the first two or three dimensions. To exclude ethnic 
    # outliers, the thresholds need to be set around the cluster of population 
    # of interest.
	  awk '{ if ($4 <0.02 && $5 >-0.025) print $1,$2 }' !{C7_eigenvec} > EUR_PCA_merge
    '''
}

// STEP C9: Extract the homogenous ethnic group of samples from user's data ----

process extract_homogenous {
  input:
  file ethnic_cluster
  file inbed_extract
  file inbim_extract
  file infam_extract

  output:
  file "C9*" into C9

  shell:
  '''
  plink --bfile sample_variant_qc \
      --keep !{ethnic_cluster} \
      --make-bed \
      --out C9
  '''
}

// STEP C10: Create covariates based on PCA ------------------------------------

process homogenous_pca {
  input:
  file C9 
  file C10_excluded_regions

  output:
  file "covar_pca.txt" into covar_homogenous_pca

  shell:
  '''
  # Pruning in user's dataset
	plink --bfile C9 \
    --exclude !{C10_excluded_regions} \
    --indep-pairwise 50 5 0.2 \
    --out indepSNPs_1k_1
	plink --bfile C9 \
    --extract indepSNPs_1k_1.prune.in \
    --make-bed \
    --out C10_indep 

  # Perform a PCA on user's data without ethnic outliers. 
	plink --bfile C10_indep \
    --pca header \
    --out C10

	# Create covariate file including the first 3 PCs
	awk '{print $1, $2, $3, $4, $5}' C10.eigenvec > covar_pca.txt
  '''
}

// Finished!
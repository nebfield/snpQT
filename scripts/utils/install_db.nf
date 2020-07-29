log.info """\
         snpQT: make your SNPs cute 
         make_ref: Process reference data downloaded by download_db.sh
         input directory : ${params.indir}
         """
         .stripIndent()

Channel.fromPath(params.indir + '/human_g1k_v37.fasta').set{ h37 }
Channel.fromPath(params.indir + '/integrated_call_samples_v3.20130502.ALL.panel').set{ in_panel } 
Channel.fromPath(params.indir + '/thousand_genomes/*.gz').set { thousand_genomes }
Channel.fromPath(params.indir + '/PCA.exclude.regions.b37.txt').set{ exclude } 

// Note: Step F1 & F2 are taken care of by download_db.sh 
// STEP F3 ---------------------------------------------------------------------
// split multi-allelics, annotate variants with a unique identifier and remove 
// duplicates

process clean_chromosomes {
  input:
  each chr from thousand_genomes
  file h37
  // each is important to recycle h37 per chromosome  

  output:
  file "*.bcf" into cleaned_chrom

  shell:
  '''
  chr_name=$(basename !{chr} .vcf.gz) # get basename of a chromosome 
  
  bcftools norm -m-any --check-ref w -f !{h37} !{chr} | \
		bcftools annotate -I '%CHROM:%POS:%REF:%ALT' | \
		bcftools norm -Ob --rm-dup both -o $chr_name.bcf
  '''
}

// Note: formerly C1 
// STEP F4: Convert reference panel to binary plink format and merge chromosomes
process make_plink {
  input:
  file chr from cleaned_chrom

  output:
  file "${chr_name}.bed" into beds
  file "${chr_name}.bed" into bims

  shell:
  '''
  chr_name=$(basename !{chr} .bcf) # get basename of a chromosome 

  plink --bcf !{chr} \
    --allow-extra-chr 0 \
    --double-id --make-bed \
    --out $chr_name
  '''
}

process merge_plink {
  input:
  file bed from beds.collect()
  file bim from bims.collect()

  output:
  file '1kG_PCA1*' into kG_PCA1

  shell:
  '''
  find . -name "*.bim" -exec basename {} .bim \\; > mergeList.txt
  plink --merge-list mergeList.txt \
    --keep-allele-order \
    --make-bed \
    --out 1kG_PCA1
  '''
}

// Note: formerly C2
// STEP F5: STEP C2: QC on 1000 Genomes data -----------------------------------

process qc_thousand_genomes {
    publishDir params.indir, mode: 'copy', overwrite: true, \
      pattern: "1kG_PCA5.bed"
    publishDir params.indir, mode: 'copy', overwrite: true, \
      pattern: "1kG_PCA5.bim"
    publishDir params.indir, mode: 'copy', overwrite: true, \
      pattern: "1kG_PCA5.fam"
      
    input:
    file kG_PCA1
    file exclude

    output:
    file "1kG_PCA6*" into thousand_genomes_qc

    shell:
    '''
    # Remove variants based on missing genotype data.
    plink --bfile 1kG_PCA1 \
      --geno 0.1 \
      --allow-no-sex \
      --make-bed \
      --out 1kG_PCA2

    # Remove individuals based on missing genotype data.
    plink --bfile 1kG_PCA2 \
      --mind 0.02 \
      --allow-no-sex \
      --make-bed \
      --out 1kG_PCA3

    # Remove variants based on missing genotype data.
    plink --bfile 1kG_PCA3 \
      --geno 0.02 \
      --allow-no-sex \
      --make-bed \
      --out 1kG_PCA4

    # Remove variants based on MAF.
    plink --bfile 1kG_PCA4 \
      --maf 0.05 \
      --allow-no-sex \
      --make-bed \
      --out 1kG_PCA5

    # Prune variants
    # Check which r2 threshold is best (Ask Andrew)
    plink --bfile 1kG_PCA5 \
      --exclude !{exclude} \
      --indep-pairwise 50 5 0.2 \
      --out indepSNPs_1k
	  plink --bfile 1kG_PCA5 \
      --extract indepSNPs_1k.prune.in \
      --make-bed \
      --out 1kG_PCA6
    '''
} 

// Step F6: Make a racefile ---------------------------------------------------
process racefile {
  publishDir params.indir, mode: 'copy', overwrite: true 

  input:
  file in_panel

  output:
  file "*.txt" into racefiles 

  shell:
  '''
  # Make 1st racefile, using the 20130502 panel using superpopulation codes 
  # (i.e., AFR,AMR,EASN,SAS and EUR)
  awk '{print$1,$1,$3}' !{in_panel}  > super_racefile_1k.txt

  # Make 2nd racefile, using the 20130502 panel using subpopulation codes 
  awk '{print$1,$1,$2}' !{in_panel} > sub_racefile_1k.txt
  '''
}


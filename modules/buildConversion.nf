// Note: Steps A1 - A3 taken care of by download_db.sh 


// Step A4: Create a dictionary file ------------------------------------------
process dictionary {
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

// Note: Step A5 is taken care of by install_db.sh
// STEP A6: Change the chr ids ------------------------------------------------
process num_to_chr {
  input:
  path(in_vcf)
  path(chr_map)
  
  output:
  path "out.vcf.gz", emit: vcf

  shell:
  '''
  bcftools annotate --rename-chr !{chr_map} !{in_vcf} \
    -Oz -o out.vcf.gz
  '''
}

// STEP A7: Run liftOver to map genome build -----------------------------------
process liftover {
  
  publishDir "${params.results}/convertBuild/files/", mode: 'copy'
  
  input:
  path(vcf)
  path(hg)
  path(chain)
  path(dict)
  
  output:
  path "out.vcf", emit: vcf

  shell:
  '''
  # !{dict} unused but needed to stage in file
  picard -Xmx16g LiftoverVcf \
    -I !{vcf} \
    -O out.vcf \
    -CHAIN !{chain} \
    -REJECT rejected_variants.vcf \
    -R !{hg}
  '''
}

// Note: A8 is combined here
// STEP A9: Reverse Chr1To1 ---------------------------------------------------
process chr_to_num {
  
  publishDir "${params.results}/convertBuild/files/", mode: 'copy'
  
  input:
  path(vcf)
  path(chr_map)
  
  output:
  path "out.vcf.gz", emit: vcf

  shell:
  '''
  awk '{print $2 "\t" $1}' !{chr_map} > Chr1To1.txt

  # Change the chromosome ids again
  bcftools annotate --rename-chr Chr1To1.txt !{vcf} -Oz -o out.vcf.gz
  '''
}

// STEP A10: Convert VCF to PLINK format ---------------------------------------
process vcf_to_plink {

  publishDir "${params.results}/convertBuild/files/", mode: 'copy'
  
  input:
  path(vcf)

  output:
  path "converted.bed", emit: bed
  path "converted.bim", emit: bim
  path "converted.fam", emit: fam

  shell:
  '''
  # Convert VCF to PLINK format
  plink2 --vcf !{vcf} \
    --max-alleles 2 \
    --chr 1-22 XY \
    --make-bed \
    --out converted 
  '''
}

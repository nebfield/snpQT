// Note: Steps A1 - A3 taken care of by download_db.sh 
// Step A4: Create a dictionary file ------------------------------------------
process dictionary {
  input:
  path(hg19)

  output:
  path "hg19.fa.gz.dict", emit: dict
  
  shell:
  '''
  # make it so!
  picard \
    CreateSequenceDictionary \
    -R !{hg19} \
    -O hg19.fa.gz.dict
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
  input:
  path(vcf)
  path(hg19)
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
    -R !{hg19}
  '''
}

// Note: A8 is combined here
// STEP A9: Reverse Chr1To1 ---------------------------------------------------
process chr_to_num {
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
  input:
  path(vcf)

  output:
  path "converted.bed", emit: bed
  path "converted.bim", emit: bim
  path "converted.fam", emit: fam

  shell:
  '''
  # Convert VCF to PLINK format
  plink --vcf !{vcf} \
    --keep-allele-order \
    --allow-extra-chr \
    --chr 1-22 X Y XY MT \
    --make-bed \
    --out converted 
  '''
}

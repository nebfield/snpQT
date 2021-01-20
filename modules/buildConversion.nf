// Note: Steps A1 - A3 taken care of by download_db.sh 
// Step A4: Create a dictionary file ------------------------------------------
process dictionary {
  input:
  path(ch_db)

  output:
  path "hg19.fa.gz.dict", emit: dict
  
  shell:
  '''
  # get specific files from db
  picard=!{ch_db}"/picard.jar"
  hg19=!{ch_db}"/hg19.fa.gz"

  echo $hg19 $picard
  
  # make it so!
  java -Dpicard.useLegacyParser=false \
    -Xmx8g \
    -jar $picard \
    CreateSequenceDictionary \
    -R $hg19 \
    -O hg19.fa.gz.dict
  '''
}

// Note: Step A5 is taken care of by install_db.sh
// STEP A6: Change the chr ids ------------------------------------------------
process num_to_chr {
  input:
  path(in_vcf)
  path(ch_db)
  
  output:
  path "out.vcf.gz", emit: vcf

  shell:
  '''
  # get specific file from db
  chr_map=!{ch_db}"/1toChr1.txt"

  bcftools annotate --rename-chr $chr_map !{in_vcf} \
    -Oz -o out.vcf.gz
  '''
}

// STEP A7: Run liftOver to map genome build -----------------------------------
process liftover {
  input:
  path(vcf)
  path(ch_db)
  path(dict)
  
  output:
  path "out.vcf", emit: vcf

  shell:
  '''
  # get specific files from db
  picard=!{ch_db}"/picard.jar"
  hg19=!{ch_db}"/hg19.fa.gz"
  chain=!{ch_db}"/hg38ToHg19.over.chain"
  cp $hg19 . # TODO: move this to the database install?
  java -Dpicard.useLegacyParser=false \
    -Xmx16g \
    -jar $picard LiftoverVcf \
    -I !{vcf} \
    -O out.vcf \
    -CHAIN $chain \
    -REJECT rejected_variants.vcf \
    -R hg19.fa.gz
  '''
}

// Note: A8 is combined here
// STEP A9: Reverse Chr1To1 ---------------------------------------------------
process chr_to_num {
  input:
  path(vcf)
  path(ch_db)
  
  output:
  path "out.vcf.gz", emit: vcf

  shell:
  '''
  # get specific files from db
  chr_map=!{ch_db}"/1toChr1.txt"
  
  awk '{print $2 "\t" $1}' $chr_map > Chr1To1.txt

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

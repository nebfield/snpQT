log.info """\
         snpQT: human genome build conversion  
         """
         .stripIndent()

params.outdir = "$baseDir/../../results"
outdir = params.outdir + '/buildConversion/'

Channel.fromPath(params.infile).set{ in_vcf }
Channel.fromPath("$SNPQT_DB_DIR/picard.jar").set{ picard_jar }
picard_jar.into{ picard_dictionary; picard_liftover}
Channel.fromPath("$SNPQT_DB_DIR/hg19.fa.gz").set{ hg19 }
hg19.into{ hg19_dictionary; hg19_liftover }
Channel.fromPath("$SNPQT_DB_DIR/hg38ToHg19.over.chain").set{ chain }
Channel.fromPath("$SNPQT_DB_DIR/1toChr1.txt").set{ chr }
chr.into { chr_map; chr_map_back }

// Note: Steps A1 - A3 taken care of by download_db.sh 
// Step A4: Create a dictionary file ------------------------------------------
process dictionary {
  container 'openjdk:8'
  // publishDir params.outdir, mode: 'copy', overwrite: true 

  input:
  file hg19_dictionary
  file picard_dictionary

  output:
  file "hg19.fa.gz.dict" into dict

  shell:
  '''
  # make it so!
  java -Dpicard.useLegacyParser=false \
    -Xmx8g \
    -jar !{picard_dictionary} \
    CreateSequenceDictionary \
    -R !{hg19_dictionary} \
    -O hg19.fa.gz.dict
  '''
}

// Note: Step A5 is taken care of by install_db.sh
// STEP A6: Change the chr ids ------------------------------------------------
process num_to_chr {
  input:
  file in_vcf
  file chr_map

  output:
  file "dataset_1.vcf.gz" into dataset_1

  shell:
  '''
  bcftools annotate --rename-chr !{chr_map} !{in_vcf} \
    -Oz -o dataset_1.vcf.gz
  '''
}

// STEP A7: Run liftOver to map genome build -----------------------------------
process liftover {
  container 'openjdk:8'

  input:
  file dataset_1
  file picard_liftover
  file chain 
  file hg19_liftover
  file dict

  output:
  file "dataset_2.vcf" into dataset_2

  shell:
  '''
  java -Dpicard.useLegacyParser=false \
    -Xmx16g \
    -jar !{picard_liftover} LiftoverVcf \
    -I dataset_1.vcf.gz \
    -O dataset_2.vcf \
    -CHAIN !{chain} \
    -REJECT rejected_variants.vcf \
    -R !{hg19_liftover}
  '''
}

// Note: A8 is combined here
// STEP A9: Reverse Chr1To1 ---------------------------------------------------
process chr_to_num {
  input:
  file dataset_2
  file chr_map_back

  output:
  file "dataset_3*" into dataset_3 

  shell:
  '''
  awk '{print $2 "\t" $1}' !{chr_map_back} > Chr1To1.txt

  # Change the chromosome ids again
  bcftools annotate --rename-chr Chr1To1.txt dataset_2.vcf -Oz -o dataset_3.vcf.gz
  '''
}

// STEP A10: Convert VCF to PLINK format ---------------------------------------
process vcf_to_plink {
  publishDir outdir, mode: 'copy', overwrite: true, pattern: "*.bed"
  publishDir outdir, mode: 'copy', overwrite: true, pattern: "*.bim"
  publishDir outdir, mode: 'copy', overwrite: true, pattern: "*.fam"

  input:
  file dataset_3

  output:
  file "dataset_4.*"

  shell:
  '''
  # Convert VCF to PLINK format
  plink --vcf dataset_3.vcf.gz \
    --keep-allele-order \
    --allow-extra-chr \
    --chr 1-22 X Y XY MT \
    --make-bed --out dataset_4 
  '''
}

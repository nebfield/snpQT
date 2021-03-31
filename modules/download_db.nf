// Note: step C2 from population stratification module

process qc_ref_data {
  input:
  path(thousand_pgen)
  path(thousand_psam)
  path(thousand_pvar)
  path(h37)
  path(exclude_region)

  output:
  path "all_phase3_1.bed", emit: bed
  path "all_phase3_1.bim", emit: bim
  path "all_phase3_1.fam", emit: fam
  path "h37_squeezed.fasta", emit: h37
  path "h37_squeezed.fasta.fai", emit: h37_idx
  
  shell:
  '''
  # -q to suppress warnings, which make nextflow unhappy
  # readlink to fix symlink gzip problem
  # || true to force exit code 0 
  gunzip --quiet -dc $(readlink !{h37}) > human_g1k_v37.fasta || true
  # fix formatting in the file or plink explodes
  tr -s ':' < human_g1k_v37.fasta | tr -s '\n' > h37_squeezed.fasta
  samtools faidx h37_squeezed.fasta
  
  mv !{thousand_psam} all_phase3.psam # fix psam name
  plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen
  plink2 --pfile 'vzs' all_phase3 --chr 1-22 XY --make-pfile 'vzs' --out all_phase3_1

  

  # Remove duplicates
  plink2 --pfile 'vzs' all_phase3_1 \
    --rm-dup force-first \
    --make-pgen 'vzs'\
    --out all_phase3_2
  # Remove multi-allelic variants
  plink2 --pfile 'vzs' all_phase3_2 \
    --max-alleles 2 \
    --make-pgen 'vzs'\
    --out all_phase3_1
  # Remove variants based on MAF
  plink2 --pfile 'vzs' all_phase3_1 \
    --maf 0.05 \
	--set-missing-var-ids @:#:\\$r:\\$a \
    --make-bed \
    --out all_phase3_2
  # Prune variants
  plink2 --bfile all_phase3_2 \
    --exclude !{exclude_region} \
    --indep-pairwise 50 5 0.2 \
    --out indepSNPs_1k_allphase
  plink2 --bfile all_phase3_2 \
    --extract indepSNPs_1k_allphase.prune.in \
    --make-bed \
    --out all_phase3_1
  '''
}

process decompress {
  input:
  path(x)

  output:
  path("${x.baseName}"), emit: file

  shell:
  '''
  bgzip -d !{x}
  '''
}

process index {
  input:
  path(x)

  output:
  path("*.tbi"), emit: idx

  shell:
  '''
  tabix !{x}
  '''
}

process qc {
  input:
  tuple val(chr), path(vcf)

  output:
  path "*.vcf.gz", emit: vcf
  
  shell:
  '''
  bcftools norm --rm-dup both  !{vcf} -Oz -o !{chr}.vcf.gz
  '''
}

process annotate_ids {
  input:
  path(vcf)

  output:
  path "*.vcf.gz", emit: vcf
  
  shell:
  '''
  bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT'  !{vcf} -Oz -o !{chr}.vcf.gz
  '''
}


process unzip_shapeit4 {
  input:
  path(x)

  output:
  path "genetic_maps.b37.tar.gz", emit: maps

  shell:
  '''
  unzip !{x}
  mv shapeit4-4.2.0/maps/genetic_maps.b37.tar.gz .
  '''
}

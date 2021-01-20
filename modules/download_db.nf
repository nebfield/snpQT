// Note: step C2 from population stratification module

process qc_ref_data {
  input:
  path(thousand_pgen)
  path(thousand_psam)
  path(thousand_pvar)
  path(h37)
  path(exclude_region)

  output:
  path "all_phase3_10.bed", emit: bed
  path "all_phase3_10.bim", emit: bim
  path "all_phase3_10.fam", emit: fam
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
  
  mv !{thousand_psam} all_phase3.psam # fix dropbox url ?dl=1
  mv !{thousand_pvar} all_phase3.pvar.zst
  plink2 --zst-decompress !{thousand_pgen} > all_phase3.pgen
  plink2 --pfile 'vzs' all_phase3 --chr 1-22 --make-pfile --out all_phase3_1

  plink2 --pfile all_phase3_1 \
    --normalize 'list' \
    --fa h37_squeezed.fasta \
    --make-pgen \
    --out all_phase3_3
  # Remove duplicates
  plink2 --pfile all_phase3_3 \
    --rm-dup force-first \
    --make-pgen \
    --out all_phase3_4
  # Remove multi-allelic variants
  plink2 --pfile all_phase3_4 \
    --max-alleles 2 \
    --make-pgen \
    --out all_phase3_5
  # Remove variants based on missing genotype data
  plink2 --pfile all_phase3_5 \
    --geno 0.1 \
    --make-pgen \
    --out all_phase3_6
  # Remove individuals based on missing genotype data.
  plink2 --pfile all_phase3_6 \
    --mind 0.02 \
    --make-pgen \
    --out all_phase3_7
  # Remove variants based on missing genotype data.
  plink2 --pfile all_phase3_7 \
    --geno 0.02 \
    --make-pgen \
    --out all_phase3_8
  # Remove variants based on MAF and prune
  plink2 --pfile all_phase3_8 \
    --maf 0.05 \
    --make-bed \
    --out all_phase3_9
  # Prune variants
  plink --bfile all_phase3_9 \
    --exclude !{exclude_region} \
    --indep-pairwise 50 5 0.2 \
    --out indepSNPs_1k_allphase
  plink --bfile all_phase3_9 \
    --extract indepSNPs_1k_allphase.prune.in \
    --make-bed \
    --out all_phase3_10
  '''
}
log.info """\
         snpQT: make your SNPs cute 
         make_ref: Process reference data downloaded by download_db.sh
         input directory : ${params.indir}
         """
         .stripIndent()

Channel.fromPath(params.indir + '/human_g1k_v37.fasta').set{ h37 }
Channel.fromPath(params.indir + '/all_phase3.pgen.zst').set { thousand_pgen }
Channel.fromPath(params.indir + '/all_phase3.pvar.zst').set { thousand_pvar }
Channel.fromPath(params.indir + '/phase3_corrected.psam').set { thousand_psam } 
Channel.fromPath(params.indir + '/PCA.exclude.regions.b37.txt').set{ exclude } 

// Note: step C2 from population stratification module

process qc {
  publishDir params.indir, mode: 'copy', overwrite: true, \
    pattern: "all_phase3_10.bed"
  publishDir params.indir, mode: 'copy', overwrite: true, \
    pattern: "all_phase3_10.fam"
  publishDir params.indir, mode: 'copy', overwrite: true, \
    pattern: "all_phase3_10.bim"

  container 'snpqt'
  
  input:
  file thousand_pgen
  file thousand_pvar
  file thousand_psam
  file h37
  file exclude

  output:
  file "all_phase3_10*"
  
  shell:
  '''
  mv !{thousand_psam} all_phase3.psam
  plink2 --zst-decompress !{thousand_pgen} > all_phase3.pgen
  plink2 --pfile 'vzs' all_phase3 --chr 1-22 --make-pfile --out all_phase3_1
  tr -s ':' !{h37} > h37_squeezed.fasta
  # fixes left-normalisation bug unexpected character
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
    --exclude !{exclude} \
    --indep-pairwise 50 5 0.2 \
    --out indepSNPs_1k_allphase
  plink --bfile all_phase3_9 \
    --extract indepSNPs_1k_allphase.prune.in \
    --make-bed \
    --out all_phase3_10
  '''
}

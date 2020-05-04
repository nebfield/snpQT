log.info """\
         snpQT: make your SNPs cute 
         make_ref: Process reference data downloaded by download_db.sh
         input directory : ${params.indir}
         """
         .stripIndent()

in_vcf = Channel.fromPath(params.indir + '/ALL.2of4intersection.20100804.genotypes.vcf.gz')
in_panel = Channel.fromPath(params.indir + '/20100804.ALL.panel')
in_race = Channel.fromPath(params.indir + '/race_1kG.txt')

process make_bed {
    input:
    file in_vcf

    output:
    file "ALL.2of4intersection*" into thousand_genomes

    """
    plink --vcf $in_vcf --make-bed \
      --out ALL.2of4intersection.20100804.genotypes 
    """
}

process make_unique_ids {
    input:
    file thousand_genomes

    output:
    file "2of4intersection.20100804.genotypes_NMIDs*" into nmids 

    shell:
    '''
    plink --bfile ALL.2of4intersection.20100804.genotypes \
      --set-missing-var-ids @:#[b37]\\$1,\\$2 --make-bed \
      --out 2of4intersection.20100804.genotypes_NMIDs
    '''
}

process qc_thousand_genomes {
    publishDir params.indir, mode: 'copy', overwrite: true, \
      pattern: "1kG_PCA5.bed"
    publishDir params.indir, mode: 'copy', overwrite: true, \
      pattern: "1kG_PCA5.bim"
    publishDir params.indir, mode: 'copy', overwrite: true, \
      pattern: "1kG_PCA5.fam"
      
    input:
    file nmids

    output:
    file "1kG_PCA5*" into thousand_genomes_qc

    """
    # Remove variants based on missing genotype data.
    plink --bfile 2of4intersection.20100804.genotypes_NMIDs --geno 0.1 \
      --allow-no-sex --make-bed --out 1kG_PCA2
    # Remove individuals based on missing genotype data.
    plink --bfile 1kG_PCA2 --mind 0.02 --allow-no-sex --make-bed --out 1kG_PCA3
    # Remove variants based on missing genotype data.
    plink --bfile 1kG_PCA3 --geno 0.02 --allow-no-sex --make-bed --out 1kG_PCA4
    # Remove variants based on MAF.
    plink --bfile 1kG_PCA4 --maf 0.05 --allow-no-sex --make-bed --out 1kG_PCA5
    """
} 

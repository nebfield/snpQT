// TODO: update params.indir 

log.info """\
         snpQT: make your SNPs cute 
         make_db: Process database downloaded by download_db.sh
         input directory : ${params.infile}
         """
         .stripIndent()

Channel
    .fromPath( params.infile )
    .ifEmpty { error "Cannot find: ${params.infile}" }
    .set { in_file } 

process make_vcf {
    echo true
    container 'snpqt'

    input:
    file in_file

    output:
    file "ALL.2of4intersection*" into thousand_genomes
    """
    plink --vcf $in_file --make-bed --out \
      ALL.2of4intersection.20100804.genotypes &>/dev/null
    """
}
process make_unique_ids {
    echo true
    container 'snpqt'

    input:
    file thousand_genomes

    output:
    file "2of4intersection.20100804.genotypes_NMIDs*" into nmids 

    """
    make_unique_ids.sh &>/dev/null
    echo 'Make unique IDs:' && \
      grep 'pass' 2of4intersection.20100804.genotypes_NMIDs.log
    """
}

process qc_thousand_genomes {
    echo true
    container 'snpqt'
    publishDir '../data/', mode: 'copy'
    input:
    file nmids

    output:
    file "1kG_MDS3*" into thousand_genomes_qc

    """
    plink -bfile 2of4intersection.20100804.genotypes_NMIDs --geno 0.02 \
      --allow-no-sex --make-bed --out 1kG_MDS1 &>/dev/null
    plink --bfile 1kG_MDS1 --mind 0.02 --allow-no-sex --make-bed \
      --out 1kG_MDS2 &>/dev/null
    plink --bfile 1kG_MDS2 --maf 0.05 --allow-no-sex --make-bed \
      --out 1kG_MDS3 &>/dev/null
    echo "QC: " && grep 'pass' 1kG_MDS3.log
    """
} 

#!/bin/bash

impute_genotypes_url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{.}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --output ALL.chr{.}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
impute_idx_url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{.}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi --output ALL.chr{.}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"

LIBRARY_DIR="$SNPQT_DB_DIR"
mkdir -p $LIBRARY_DIR

cd $LIBRARY_DIR

seq 1 22 > chrom_list.txt

if [ ! -e "impute.download.complete" ]
then
    parallel -a chrom_list.txt -j 6 --progress --resume-failed --joblog log curl $impute_genotypes_url 	
    parallel -a chrom_list.txt -j 6 --progress --resume-failed --joblog log curl $impute_idx_url
    touch impute.download.complete
fi

# note: human_g1k_v37.fasta.gz taken care of by download_db.sh


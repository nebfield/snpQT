#!/usr/bin/env bash

# Download refernce data

set -u  
set -e  

# URL list ---------------------------------------------------------------------

# build conversion
picard_url="https://github.com/broadinstitute/picard/releases/download/2.22.4/picard.jar"
chain_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
hg19_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"

# reference genome
human_url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"
human_fai_url="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai"

# thousand genome
panel_url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
thousand_vcf_url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{.}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
thousand_tabix_url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{.}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"

LIBRARY_DIR="$SNPQT_DB_DIR"
mkdir -p $LIBRARY_DIR

cd $LIBRARY_DIR

if [ ! -e "download.complete" ]
    then
    # small files --------------------------------------------------------------
    curl -s -O $panel_url
    curl -s -O -L $picard_url
    curl -s -O $chain_url
    curl -s -O $human_fai_url 
    bgzip -d hg38ToHg19.over.chain.gz 
    # bigger files -------------------------------------------------------------
    echo -n "Downloading GRCh37 (~1GB)..."
    curl -s -O $human_url
    echo " finished."
    bgzip -d human_g1k_v37.fasta.gz
    echo -n "Downloading hg19 (~1GB)..."
    curl -s -O $hg19_url
    echo " finished."
    bgzip -d hg19.fa.gz
    touch "download.complete"
fi

if [ ! -e "big.download.complete" ]
    then
    echo -n "Downloading thousand genomes genotypes..."
    seq 1 22 > chrom_list.txt
    mkdir -p thousand_genomes/ && cd thousand_genomes 
    parallel -a ../chrom_list.txt --progress --resume-failed --joblog vcf_log \
      curl -s -O $thousand_vcf_url
    parallel -a ../chrom_list.txt --progress --resume-failed --joblog tabix_log \
      curl -s -O $thousand_tabix_url 
    rm vcf_log tabix_log
    echo " finished."
    cd .. 
    touch "big.download.complete" 
fi

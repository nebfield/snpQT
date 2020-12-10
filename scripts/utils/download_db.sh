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

# 1K genome data from plink (pop strat)
psam_url="https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam?dl=1"
pgen_url="https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1"
pvar_url="https://www.dropbox.com/s/op9osq6luy3pjg8/all_phase3.pvar.zst?dl=1"

# dbSNP (step D8)
dbsnp_url="ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz"
dbsnp_index_url="ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz.tbi"

# shapeit4 
shapeit4_maps_url="https://github.com/odelaneau/shapeit4/blob/master/maps/genetic_maps.b37.tar.gz?raw=true"

LIBRARY_DIR="$SNPQT_DB_DIR"
mkdir -p $LIBRARY_DIR

cd $LIBRARY_DIR

if [ ! -e "download.complete" ]
then
    # small files --------------------------------------------------------------
    curl -s -O -L $picard_url
    curl -s -O $chain_url
    curl -s -O $human_fai_url 
    bgzip -d hg38ToHg19.over.chain.gz 
    curl -sL $shapeit4_maps_url -o genetic_maps.b37.tar.gz
    # bigger files -------------------------------------------------------------
    echo -n "Downloading GRCh37 (~1GB)..."
    curl -s -O $human_url
    echo " finished."
    bgzip -d human_g1k_v37.fasta.gz
    echo -n "Downloading hg19 (~1GB)..."
    curl -s -O $hg19_url
    echo " finished."
    touch "download.complete"
fi

if [ ! -e "big.download.complete" ]
then
    echo -n "Downloading thousand genomes genotypes..."
    curl -Lso phase3_corrected.psam $psam_url
    curl -Lso all_phase3.pgen.zst $pgen_url
    curl -Lso all_phase3.pvar.zst $pvar_url
    echo " finished."
    echo -n "Downloading dbSNP (~15GB)..."
    curl -s -O $dbsnp_url
    curl -s -O $dbsnp_index_url
    echo " finished."
    touch "big.download.complete" 
fi

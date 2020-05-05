#!/usr/bin/env bash

# Download 1000 genomes data for use with snpQT population stratification

set -u  
set -e  

LIBRARY_DIR="$SNPQT_DB_DIR"
gt_url="ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz"
panel_url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/20100804.ALL.panel"
human_url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"

# URL list ---------------------------------------------------------------------
# vcf sanity checking

picard_url="https://github.com/broadinstitute/picard/releases/download/2.22.4/picard.jar"
chain_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
hg19_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"

mkdir -p $LIBRARY_DIR
cd $LIBRARY_DIR

if [ ! -e "small.download.complete" ]
    then
    wget $panel_url
    wget $picard_url
    wget $chain_url
    touch "small.download.complete"
fi

if [ ! -e "medium.download.complete" ]
    then
    echo "Downloading GRCh37 (~1GB)..."
    wget $human_url
    echo " finished."
    echo "Downloading hg19 (~1GB)..."
    wget $hg19_url
    echo " finished."
    touch "medium.download.complete"
fi

if [ ! -e "big.download.complete" ]
    then
    echo "Downloading thousand genomes genotypes (~60GB)..."
    wget $gt_url 
    echo " finished."
    touch "big.download.complete"
fi

# decompress downloaded files
bgzip -d human_g1k_v37.fasta.gz
bgzip -d hg19.fa.gz
bgzip -d hg38ToHg19.over.chain.gz

# index 
samtools faidx hg19.fa

# run picard 
java -Dpicard.useLegacyParser=false -jar picard.jar \
      CreateSequenceDictionary \
      -R hg19.fa \
      -O hg19.fa.dict

# Convert population codes into superpopulation codes (i.e., AFR,AMR,ASN,
# and EUR).
awk '{print$1,$1,$2}' 20100804.ALL.panel > 1kG_race.txt
sed -i 's/JPT/ASN/g' 1kG_race.txt
sed -i 's/ASW/AFR/g' 1kG_race.txt
sed -i 's/CEU/EUR/g' 1kG_race.txt
sed -i 's/CHB/ASN/g' 1kG_race.txt
sed -i 's/CHD/ASN/g' 1kG_race.txt
sed -i 's/YRI/AFR/g' 1kG_race.txt
sed -i 's/LWK/AFR/g' 1kG_race.txt
sed -i 's/TSI/EUR/g' 1kG_race.txt
sed -i 's/MXL/AMR/g' 1kG_race.txt
sed -i 's/GBR/EUR/g' 1kG_race.txt
sed -i 's/FIN/EUR/g' 1kG_race.txt
sed -i 's/CHS/ASN/g' 1kG_race.txt
sed -i 's/PUR/AMR/g' 1kG_race.txt

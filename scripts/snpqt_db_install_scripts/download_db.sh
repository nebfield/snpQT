#!/usr/bin/env bash

# Download 1000 genomes data for use with snpQT population stratification

set -u  
set -e  

LIBRARY_DIR="$SNPQT_DB_DIR/thousand_genomes"
gt_url="ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz"
panel_url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/20100804.ALL.panel"

mkdir -p $LIBRARY_DIR
cd $LIBRARY_DIR

if [ ! -e "download.complete" ]
    then
    echo "Downloading thousand genomes population information..."
    wget $panel_url
    echo " finished."
    echo "Downloading thousand genomes genotypes (~60GB)..."
    wget $gt_url 
    echo " finished."
    touch "download.complete"
fi

# Convert population codes into superpopulation codes (i.e., AFR,AMR,ASN,
# and EUR).
awk '{print$1,$1,$2}' 20100804.ALL.panel > race_1kG.txt
sed -i 's/JPT/ASN/g' race_1kG.txt
sed -i 's/ASW/AFR/g' race_1kG.txt
sed -i 's/CEU/EUR/g' race_1kG.txt
sed -i 's/CHB/ASN/g' race_1kG.txt
sed -i 's/CHD/ASN/g' race_1kG.txt
sed -i 's/YRI/AFR/g' race_1kG.txt
sed -i 's/LWK/AFR/g' race_1kG.txt
sed -i 's/TSI/EUR/g' race_1kG.txt
sed -i 's/MXL/AMR/g' race_1kG.txt
sed -i 's/GBR/EUR/g' race_1kG.txt
sed -i 's/FIN/EUR/g' race_1kG.txt
sed -i 's/CHS/ASN/g' race_1kG.txt
sed -i 's/PUR/AMR/g' race_1kG.txt
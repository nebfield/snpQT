#!/usr/bin/env bash

# Download the file with population information of the 1000 genomes
# dataset.
curl -O ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/20100804.ALL.panel
# The file 20100804.ALL.panel contains population codes of the individuals
# of 1000 genomes.

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


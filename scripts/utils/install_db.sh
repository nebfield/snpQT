#!/usr/bin/env bash

# QC downloaded 1000 genomes data for use with snpQT population stratification

# set -e
# set -u 

export NXF_WORK=$(realpath "../../work") # keep work directories at top level 
SNPQT_CONFIG=$(realpath "../nextflow.config")
PIPELINE_PATH=$(realpath "install_db.nf")
# IN_VCF=$SNPQT_DB_DIR'/ALL.2of4intersection.20100804.genotypes.vcf.gz'
# # VCF_PATH=$(realpath $IN_VCF)

# if [ ! -f "$VCF_PATH" ]; then
#     echo "$VCF_PATH doesn't exist. Did you run download_db.sh? Is \$SNPQT_DB_DIR set?"
#     exit 1
# fi

cp $(realpath "../../data/PCA.exclude.regions.b37.txt") $SNPQT_DB_DIR/

nextflow run $PIPELINE_PATH \
  -c $SNPQT_CONFIG \
  --indir $SNPQT_DB_DIR \
  -resume


# cp $(realpath "../../db/1toChr1.txt") $SNPQT_DB_DIR/

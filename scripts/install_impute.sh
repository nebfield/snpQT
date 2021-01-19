#!/usr/bin/env bash

# QC downloaded 1000 genomes data for use with imputation
# set -e
# set -u

export NXF_WORK=$(realpath "../../work") # keep work directories at top level
SNPQT_CONFIG=$(realpath "../nextflow.config")
PIPELINE_PATH=$(realpath "install_impute.nf")

nextflow run $PIPELINE_PATH \
  -c $SNPQT_CONFIG \
  --indir $SNPQT_DB_DIR \
  -resume

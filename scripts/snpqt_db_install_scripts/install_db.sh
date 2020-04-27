#!/usr/bin/env bash

# QC downloaded 1000 genomes data for use with snpQT population stratification

set -e
set -u 

NXF_WORK=$(realpath "../../work") # keep work directories at top level 
SNPQT_CONFIG=$(realpath "../nextflow.config")

nextflow run make_ref.nf \
  -c $SNPQT_CONFIG \
  --indir $SNPQT_DB_DIR
  -resume 
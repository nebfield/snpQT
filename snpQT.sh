#!/usr/bin/env bash

# main script
# will probably eventually be in python 

# step1="01.sample_variant.nf"

# nextflow run "$step1" #--infile "../data/als.vcf.gz" --famfile "../data/plinkForDbgap12319.fam" --outdir "results" 

# Step 1 
nextflow run scripts/01_samplevariant/01.sample_variant.nf \
  --infile '../data/als_sub.vcf.gz' \
  --famfile '../data/subset.fam' \
  --highldregion '../data/highldregion_37.txt' \
  --outdir 'results'
#!/usr/bin/env bash
set -e

# main script
# will probably eventually be in python  --------------------------------------

# POSIX method of getting path of executing script
# https://stackoverflow.com/questions/630372/determine-the-path-of-the-executing-bash-script
prg=$0
if [ ! -e "$prg" ]; then
  case $prg in
    (*/*) exit 1;;
    (*) prg=$(command -v -- "$prg") || exit;;
  esac
fi
dir=$(
  cd -P -- "$(dirname -- "$prg")" && pwd -P
) || exit

# config variables
export NXF_WORK=$(realpath work) # keep work directories at top level 
SNPQT_CONFIG=$(realpath scripts/nextflow.config)

# Step 0: ignore this for now -------------------------------------------------
nextflow run scripts/00_vcf/00.vcf.nf \
    -c $SNPQT_CONFIG \
    --infile '../data/als_sub.vcf.gz' \
    --outdir "$PWD/results" \
    -resume
 
# Step 1 ----------------------------------------------------------------------
# Sample-variant QC 
nextflow run scripts/01_samplevariant/01.sample_variant.nf \
   -c $SNPQT_CONFIG \
   --infile $(realpath '../data/als_sub.vcf.gz') \
   --inbim $(realpath 'results/vcf_sanity/dataset_4.bim') \
   --inbed $(realpath 'results/vcf_sanity/dataset_4.bed') \
   --infam $(realpath '../data/subset.fam') \
   --outdir $(realpath 'results/') \
   -resume 

# Step 2 ----------------------------------------------------------------------
# Population stratification
nextflow run scripts/02_popstrat/02.pop_strat.nf \
  -c $SNPQT_CONFIG \
  --inbed $(realpath "results/sample_qc/sample_variant_qc.bed") \
  --inbim $(realpath "results/sample_qc/sample_variant_qc.bim") \
  --infam $(realpath "results/sample_qc/sample_variant_qc.fam") \
  --outdir $(realpath 'results/') \
  -resume 
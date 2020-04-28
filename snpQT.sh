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
# nextflow run scripts/00_vcf/00.vcf.nf \
#   -c $SNPQT_CONFIG \
#   --infile '../data/als_sub.vcf.gz' \
#   --outdir "$PWD/results" \
#   --ref_fasta "../data/human_g1k_v37.fasta" 

# Step 1 ----------------------------------------------------------------------
nextflow run scripts/01_samplevariant/01.sample_variant.nf \
  -c $SNPQT_CONFIG \
  --infile $(realpath '../data/als_sub.vcf.gz') \
  --famfile $(realpath '../data/subset.fam') \
  --outdir $(realpath 'results/') \
  -resume 
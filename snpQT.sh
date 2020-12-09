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

# Step 1 ----------------------------------------------------------------------
# VCF sanity checking
nextflow run scripts/01.buildConversion/01.buildConversion.nf \
    -c $SNPQT_CONFIG \
    --infile '../data/als_sub.vcf.gz' \
    --outdir "$PWD/results" \
    -resume \
    -with-report reports/buildConversion.html
 
# Step 1 ----------------------------------------------------------------------
# Sample-variant QC 
nextflow run scripts/02.mainQC/02.mainQC.nf \
   -c $SNPQT_CONFIG \
   --infile $(realpath '../data/als_sub.vcf.gz') \
   --inbim $(realpath 'results/buildConversion/dataset_4.bim') \
   --inbed $(realpath 'results/buildConversion/dataset_4.bed') \
   --infam $(realpath '../data/subset.fam') \
   --outdir $(realpath 'results/') \
   -resume \
   -with-report reports/mainQC.html

# Step 2 ----------------------------------------------------------------------
# Population stratification
nextflow run scripts/02.popStrat/02.pop_strat.nf \
  -c $SNPQT_CONFIG \
  --inbed $(realpath "results/sample_qc/sample_variant_qc.bed") \
  --inbim $(realpath "results/sample_qc/sample_variant_qc.bim") \
  --infam $(realpath "results/sample_qc/sample_variant_qc.fam") \
  --outdir $(realpath 'results/') \
  -resume \
  -with-report reports/popStrat.html

# Step 3 ----------------------------------------------------------------------
# Imputation

nextflow run scripts/03.imputation/03.imputation.nf \
  -c $SNPQT_CONFIG \
  --inbed $(realpath "results/sample_qc/sample_variant_qc.bed") \
  --inbim $(realpath "results/sample_qc/sample_variant_qc.bim") \
  --infam $(realpath "results/sample_qc/sample_variant_qc.fam") \
  --outdir $(realpath 'results/') \
  -resume \
  -with-report reports/imputation.html

# Step 4 -----------------------------------------------------------------------
# Post-imputation QC

nextflow run scripts/04.postImputation/04.postImputation.nf \
	 -c $SNPQT_CONFIG \
	 --inimp $(realpath "results/imputation/imputed_chr*") \
	 -resume \
	 with-report reports/postimputation.html
	 

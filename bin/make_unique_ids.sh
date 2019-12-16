#!/usr/bin/env bash

# Assign unique indentifiers to the SNPs with a missing rs-identifier \
# (i.e., the SNPs with ".")
plink --bfile ALL.2of4intersection.20100804.genotypes \
  --set-missing-var-ids @:#[b37]$1,$2 --make-bed \
  --out 2of4intersection.20100804.genotypes_NMIDs


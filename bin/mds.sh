#!/usr/bin/env bash

# Perform MDS on plink data anchored by 1000 Genomes data.
# Using a set of pruned SNPs

plink --bfile MDS_merge --extract independent_SNPs.prune.in --genome \
	--out MDS_merge &>/dev/null
plink --bfile MDS_merge --read-genome MDS_merge.genome --cluster \
	--mds-plot 10 --out MDS_merge &>/dev/null

awk '{ if ($4 <0.02 && $5 >-0.025) print $1,$2 }' MDS_merge.mds > EUR_MDS_merge

plink --bfile missing --keep EUR_MDS_merge --make-bed --out plink_6

# Create covariates based on MDS.
# Perform an MDS ONLY on HapMap data without ethnic outliers. The values of the
# 10 MDS dimensions are subsequently used as covariates in the association
# analysis in the third tutorial.

plink --bfile plink_6 --extract independent_SNPs.prune.in --genome \
	--out plink_6
plink --bfile plink_6 --read-genome plink_6.genome --cluster \
	--mds-plot 10 --out plink_6_mds

# Change the format of the .mds file into a plink covariate file.

awk '{print$1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' \
	plink_6_mds.mds > covar_mds.txt



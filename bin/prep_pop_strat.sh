#!/usr/bin/env bash

# Extract the variants present in HapMap dataset from the 1000 genomes dataset.
awk '{print$2}' missing.bim > plink_5_SNPs.txt plink --bfile 1kG_MDS3 \
	--extract plink_5_SNPs.txt --make-bed --out 1kG_MDS3_plink_5 &>/dev/null

# Extract the variants present in 1000 Genomes dataset from the HapMap dataset.
awk '{print$2}' 1kG_MDS3_plink_5.bim > 1kG_MDS3_SNPs.txt plink \
	--bfile missing --extract 1kG_MDS3_SNPs.txt --recode \
	--make-bed --out plink_MDS &>/dev/null

## The datasets must have the same build. Change the build 1000 Genomes data
build.  awk '{print$2,$4}' plink_MDS.map > build_map.txt
# build_map.txt contains one SNP-id and physical position per line.

plink --bfile 1kG_MDS3_plink_5 --update-map build_map.txt --make-bed \
	--out 1kG_MDS4_plink_5 &>/dev/null
# 1kG_MDS4_plink_5 and plink_5 now have the same build.

## Merge the HapMap and 1000 Genomes data sets

# Prior to merging 1000 Genomes data with the HapMap data we want to make sure
# that the files are mergeable, for this we conduct 3 steps: 1) Make sure the
# reference genome is similar in the HapMap and the 1000 Genomes Project
# datasets.  2) Resolve strand issues.  3) Remove the SNPs which after the
# previoItaly two steps still differ between datasets.

# The following steps are maybe quite technical in terms of commands, but we
# jItalyt compare the two data sets and make sure they correspond.

# 1) set reference genome 

awk '{print$2,$5}' 1kG_MDS4_plink_5.bim > 1kg_ref-list.txt 
plink --bfile plink_MDS --reference-allele 1kg_ref-list.txt --make-bed \
	--out plink_MDS-adj &>/dev/null

# The 1kG_MDS4_plink_5 and the plink_MDS-adj have the same reference genome for
# all SNPs.  This command will generate some warnings for impossible A1 allele
# assignment.

# 2) Resolve strand issues.  Check for potential strand issues.

awk '{print$2,$5,$6}' 1kG_MDS4_plink_5.bim > 1kGMDS4_tmp 
awk '{print$2,$5,$6}' plink_MDS-adj.bim > plink_MDS-adj_tmp 
sort 1kGMDS4_tmp plink_MDS-adj_tmp | uniq -u > all_differences.txt

## Flip SNPs for resolving strand issues.
# Print SNP-identifier and remove duplicates.

awk '{print$1}' all_differences.txt | sort -u > flip_list.txt

# Generates a file of SNPs that are the non-corresponding between the two
# files. 

# Flip the non-corresponding SNPs. 
plink --bfile plink_MDS-adj --flip flip_list.txt \
	--reference-allele 1kg_ref-list.txt --make-bed \
	--out corrected &>/dev/null

# Check for SNPs which are still problematic after they have been flipped.
awk '{print$2,$5,$6}' corrected.bim > corrected_tmp 
sort 1kGMDS4_tmp corrected_tmp | uniq -u  > uncorresponding_SNPs.txt

# This file demonstrates the differences between the files.

# 3) Remove problematic SNPs from plink file and 1000 Genomes.
awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exclusion.txt
# The command above generates a list of the SNPs which caused the differences
# between the plink file and the 1000 Genomes data sets after flipping and
# setting of the reference genome.

# Remove the problematic SNPs from both datasets.
plink --bfile corrected --exclude SNPs_for_exclusion.txt --make-bed \
	--out plink_MDS2 &>/dev/null
plink --bfile 1kG_MDS4_plink_5 --exclude SNPs_for_exclusion.txt \
	--make-bed --out 1kG_MDS5_plink_5 &>/dev/null

# Merge plink with 1000 Genomes Data.
plink --bfile plink_MDS2 --bmerge 1kG_MDS5_plink_5.bed 1kG_MDS5_plink_5.bim \
	1kG_MDS5_plink_5.fam --allow-no-sex --make-bed \
	--out MDS_merge &>/dev/null

# racefiles -------------------------------------------------------------------
# Create a racefile of your own data.
awk '{print$1,$2,"OWN"}' plink_MDS2.fam>racefile_own.txt

# Concatenate racefiles.
cat race_1kG.txt racefile_own.txt | sed -e '1i\FID IID race' > racefile.txt


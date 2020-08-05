
// Pre-imputation 
// =============================================================================

// STEP D1: Check for strand issues, ambiguous and duplicated SNPs using snpflip

// STEP D2: Remove ambiguous SNPs ---------------------------------------------

// STEP D3: Remove one of each pair of duplicated SNPs ------------------------

// STEP D4: Flip all SNPs that are on the reverse strand ----------------------

// STEP D5: Convert Plink file into VCF ---------------------------------------

// STEP D6: bgzip VCF and then VCF needs to be indexed/sorted -----------------

// STEP D7: Convert .vcf.gz file to .bcf file ---------------------------------

// STEP D8: Check and fix the REF allele --------------------------------------

// STEP D9: Sort the BCF ------------------------------------------------------

// STEP D10: Convert .bcf file to .vcf.gz file --------------------------------

// STEP D11: Index the vcf.gz -------------------------------------------------

// STEP D12: Split vcf.gz file in chromosomes ---------------------------------

// STEP D13: Index all chroms .vcf.gz -----------------------------------------

// STEP D14: Perform phasing using shapeit4 -----------------------------------

// STEP D15: Index phased chromosomes -----------------------------------------

// Imputation 
// =============================================================================

// Note: STEP D16 is taken care of by Dockerfile 
// STEP D17: Convert vcf reference genome into a .imp5 format for each chromosome

// STEP D18: Perform imputation using impute5 ---------------------------------

// Finished!
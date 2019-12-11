# snpQT: make your SNPs cute 

## General notes

* User needs to provide a fam file for sexcheck to work  
* docker requires root

## Problems 

* Converting Christina's PLINK data to VCF caused some problems with underscores in sample names
  * `tr '_' '-' < plinkForDbgap12319.fam > new.fam` # update fam file
  *  `plink -bfile plinkForDbgap12319 --recode vcf bgz --out als` # make a new bgzipped vcf
* pihat file empty? 
* High LD regions for 38 genome build? Crossmap them? 


## TODO

* use basenames for plink output 

# snpQT: make your SNPs cute 

## General notes

* User needs to provide a fam file for sexcheck to work  
* docker requires root

# Set up 

## Download and install reference databases 

* Set `SNPQT_DB_DIR`e.g. /data/projects/qc_pipeline/db
* Run `scripts/snpqt_db_install_scripts/download_db.sh` (~60GB)
* Run `scripts/snpqt_db_install_scripts/install_db.sh`

## Problems 

* Converting Christina's PLINK data to VCF caused some problems with underscores in sample names
  * `tr '_' '-' < plinkForDbgap12319.fam > new.fam` # update fam file
  *  `plink -bfile plinkForDbgap12319 --recode vcf bgz --out als` # make a new bgzipped vcf
* pihat file empty? 

## Christina's Updates 
* Remove the --remove HighLDRegions.txt flag from the pruning step. So, we do not need a high ld regions file for each build. 
* Changed the --indep-pairwise thresholds from 50 5 0.2 to 50 5 0.5 at the same step/command
* Clone my repository where I keep the QC generalised pipeline so that you see all the changes I am making instantly



## TODO

* use basenames for plink output 

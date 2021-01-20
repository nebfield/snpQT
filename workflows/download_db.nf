// Database download workflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {qc_ref_data} from '../modules/download_db.nf' 
include {decompress} from '../modules/download_db.nf'
include {index} from '../modules/download_db.nf'

// workflow component for snpqt pipeline
workflow download_core {
  main:
    // nextflow paths can use https & ftp?! neat!
    // https paths require ? to be escaped with \\
    Channel
      .fromPath("https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam\\?dl=1", checkIfExists: true)
      .set{thousand_psam}
    Channel
      .fromPath("https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst\\?dl=1", checkIfExists: true)
      .set{thousand_pgen}
    Channel
      .fromPath("https://www.dropbox.com/s/op9osq6luy3pjg8/all_phase3.pvar.zst\\?dl=1", checkIfExists: true)
      .set{thousand_pvar}
    Channel
      .fromPath("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz", checkIfExists: true)
      .set{h37}
    Channel
      .fromPath("$baseDir/bootstrap/PCA.exclude.regions.b37.txt")
      .set{exclude_regions}
      
    // processing of core data 
    qc_ref_data(thousand_pgen, thousand_psam, thousand_pvar, h37, exclude_regions)

    // misc files
    // build conversion
    Channel
      .fromPath("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz", checkIfExists: true)
      .set{hg19}
    Channel
      .fromPath("$baseDir/bootstrap/1toChr1.txt", checkIfExists: true)
      .set{chr}
    Channel
      .fromPath("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz", checkIfExists: true)
      .set{chain}
    decompress(chain)
    // population stratification
    Channel
      .fromPath("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel", checkIfExists: true)
      .set{panel}
      
    // publish to db directory    
    println "Downloading core database files..."
    qc_ref_data.out.bed
      .concat(qc_ref_data.out.bim, qc_ref_data.out.fam, qc_ref_data.out.h37, qc_ref_data.out.h37_idx, exclude_regions, hg19, chr, decompress.out.file, panel)
      .collectFile(storeDir: "$baseDir/db/")
}

workflow download_impute {
  main:
    // nextflow paths
    Channel
      .fromPath("https://github.com/odelaneau/shapeit4/blob/master/maps/genetic_maps.b37.tar.gz\\?raw=true", checkIfExists: true)
      .set{shapeit4_maps}
    Channel
      .fromPath("ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz", checkIfExists: true)
      .set{dbsnp}
    Channel
      .fromPath("ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz.tbi", checkIfExists: true)
      .set{dbsnp_idx}

    // TODO: figure out a way to nicely download these files & idxs
    Channel
      .of(1..22)
      .map{ chrom -> "'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chrom.toString()}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'" }

    // TODO: horrible method
    Channel
      .fromPath(['ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz','ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'])
      .set{ urls }
    index(urls)
    // publish to db directory
    println "Downloading database files for imputation, this might take a while! Go and have a cup of tea :)"
    shapeit4_maps
      .concat(dbsnp, dbsnp_idx, urls, index.out.idx.collect())
      .collectFile(storeDir: "$baseDir/db/impute")
}
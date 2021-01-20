// Database download workflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {qc_ref_data} from '../modules/download_db.nf' 

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
    // population stratification
    Channel
      .fromPath("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel", checkIfExists: true)
      .set{panel}
      
    // publish to db directory
    qc_ref_data.out.bed
      .concat(qc_ref_data.out.bim, qc_ref_data.out.fam, qc_ref_data.out.h37, qc_ref_data.out.h37_idx, exclude_regions, hg19, chr, chain, panel)
      .collectFile(storeDir: "$baseDir/db/")
}

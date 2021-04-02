// Database download workflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {qc_ref_data} from '../modules/download_db.nf' // A1
include {decompress; decompress as decompress_otherBuild} from '../modules/download_db.nf' // A2
include {qc} from '../modules/download_db.nf' // A3
include {unzip_shapeit4} from '../modules/download_db.nf' // A4
include {annotate_ids} from '../modules/download_db.nf' // A5

// workflow component for snpqt pipeline
workflow download_core {
  main:
    // nextflow paths can use https & ftp?! neat!
    // https paths require ? to be escaped with \\
    Channel
      .fromPath("ftp://parrot.genomics.cn/gigadb/pub/10.5524/100001_101000/100516/phase3_corrected.psam", checkIfExists: true)
      .set{thousand_psam}
    Channel
      .fromPath("ftp://parrot.genomics.cn/gigadb/pub/10.5524/100001_101000/100516/all_phase3.pgen.zst", checkIfExists: true)
      .set{thousand_pgen}
    Channel
      .fromPath("ftp://parrot.genomics.cn/gigadb/pub/10.5524/100001_101000/100516/all_phase3.pvar.zst", checkIfExists: true)
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
	// from b38 to h19
    Channel
      .fromPath("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz", checkIfExists: true)
      .set{hg19}
    Channel
      .fromPath("$baseDir/bootstrap/1toChr1.txt", checkIfExists: true)
      .set{chr}
    Channel
      .fromPath("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz", checkIfExists: true)
      .set{chain_hg38to19}
    decompress(chain_hg38to19)
	// from h19 to b38
	Channel
      .fromPath("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz", checkIfExists: true)
      .set{hg38}
	Channel
      .fromPath("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", checkIfExists: true)
      .set{chain_hg19to38}
    decompress_otherBuild(chain_hg19to38)
    // population stratification
    Channel
      .fromPath("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel", checkIfExists: true)
      .set{panel}
	// pre-imputation
	Channel
      .fromPath("ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz", checkIfExists: true)
      .set{dbsnp}
    Channel
      .fromPath("ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz.tbi", checkIfExists: true)
      .set{dbsnp_idx} 
    // publish to db directory    
    println "Downloading core database files..."
    qc_ref_data.out.bed
      .concat(qc_ref_data.out.bim, qc_ref_data.out.fam, qc_ref_data.out.h37, qc_ref_data.out.h37_idx, exclude_regions, hg19, hg38, chr, decompress.out.file, decompress_otherBuild.out.file, panel, dbsnp, dbsnp_idx)
      .collectFile(storeDir: "$baseDir/db/")
}

workflow download_impute {
  main:
    // nextflow paths
    Channel
      .fromPath("https://github.com/odelaneau/shapeit4/archive/v4.2.0.zip", checkIfExists: true)
      .set{shapeit4_maps}
    unzip_shapeit4(shapeit4_maps)

    Channel
      .of(1..22)
      .map{ chrom -> "'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chrom.toString()}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz'" }

    Channel
      .fromPath(['ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz','ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz'])
      .map{ f -> [f.baseName.find(/\b(chr[a-zA-Z0-9]+)/), f] } // word boundary e.g. chr1, chrX
      .set{ urls }

    Channel
      .fromPath(['ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi','ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi', 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz.tbi'])
      .map{ f -> [f.baseName.find(/\b(chr[a-zA-Z0-9]+)/) + '.vcf.gz.tbi', f] } // word boundary e.g. chr1, chrX
      .set{ idx_urls }
      // we don't actually join indexes by chr so don't need chr to match
      // need .vcf.gz.tbi extension for .collectFile()

    Channel
      .fromPath("$baseDir/db/h37_squeezed.fasta", checkIfExists: true)
      .set{g37}

    qc(urls)
    annotate_ids(qc.out.vcf)
	
    // publish to db directory
    println "Downloading database files for imputation, this might take a while! Go and have a cup of tea :)"
    unzip_shapeit4.out.maps
      .concat(annotate_ids.out.vcf, idx_urls)
      .collectFile(storeDir: "$baseDir/db/impute")
}

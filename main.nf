#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {printHelp} from './modules/help.nf'

// import subworkflows
include {buildConversion} from './workflows/buildConversion.nf'
include {sample_qc} from './workflows/sample_qc.nf'
include {variant_qc} from './workflows/variant_qc.nf'
include {popStrat} from './workflows/popStrat.nf'
include {imputation} from './workflows/imputation.nf'
include {postImputation} from './workflows/postImputation.nf'
include {gwas} from './workflows/gwas.nf'
include {download_core} from './workflows/download_db.nf'
include {download_impute} from './workflows/download_db.nf'

if (params.help) {
  printHelp()
  System.exit(0)
}

println """
=================================================================
snpQT is ready to make your single-nucleotide polymorphisms cute!
v1.0, January 2020
Parameters in effect:
${params}
=================================================================
        """.stripIndent()

// throw errors on invalid workflow combinations --------------------------
if (params.convertBuild) {
  if (!params.vcf) {
    println("Please supply a vcf.gz file for build conversion with --vcf")
    println("Use --help to print help")
    System.exit(1)
  }
  if (!params.fam) {
    println("Please supply a .fam file for build conversion with --fam")
    println("Use --help to print help")
    System.exit(1)
  }
} else if (!params.convertBuild && !params.download_db) {
  if (params.vcf) {
    println("--vcf only compatible with --convertBuild")
    println("Please supply plink input files with --bed --bim --fam")
    println("Use --help to print help")
    System.exit(1)
  }
  if (!params.fam) {
    println("Missing --fam input")
    println("Use --help to print help")
    System.exit(1)
    }
}

if (params.qc) {
 if (params.bed && !params.bim || !params.bed && params.bim ) {
    println("--bed and --bim must be supplied together")
    println("Use --help to print help")
    System.exit(1)
  }
} else if (!params.qc && params.popStrat) {
  println("--popStrat requires --qc")
  println("Use --help to print help")
  System.exit(1)
} else if (!params.qc && params.impute) {
  println("--impute requires --qc")
  println("Use --help to print help")
  System.exit(1)
}

if (params.gwas) {
  if (!params.qc) {
    println("GWAS module requires qc")
    println("Please rerun with --qc")
    println("Use --help to print help")
    System.exit(1)
  }
}

// throw errors on dumb mistakes I've made
if (params.convertBuild || params.qc ) {
  if (file(params.fam).getExtension() != "fam") {
    println("Your fam file doesn't have a .fam extension. Are you sure about that?")
    println("Please rename your fam file")
    System.exit(1)
  }
}

if (params.popstrat) {
  println("I think you mean --popStrat not --popstrat. Please try again with --popStrat")
  System.exit(1)
}

// main workflow ----------------------------------------------------------
workflow {
  // set up input channels
  if ( params.convertBuild ) {
    Channel
      .fromPath(params.vcf, checkIfExists: true)
      .set{ ch_vcf }
  }
  
  if ( !params.convertBuild && params.qc ) {
    Channel
      .fromPath(params.bed, checkIfExists: true)
      .set{ ch_bed }
    Channel
      .fromPath(params.bim, checkIfExists: true)
      .set{ ch_bim }
    Channel
      .fromPath(params.fam, checkIfExists: true)
      .set{ ch_fam }
  } else if (params.convertBuild && params.qc) {
      Channel
        .fromPath(params.fam, checkIfExists: true)
        .set{ ch_fam }
  }

 if (!params.popStrat && params.gwas ) {
  Channel.fromPath("$baseDir/bootstrap/covar.txt").set{ dummy_covar }
 }

  main:
    if ( params.download_db == "core" ) {
      download_core()
    } else if (params.download_db == "impute") {
      download_impute()
    } 

    // very messy, but it works!
    // workflow with build conversion
    if ( params.convertBuild) {
      buildConversion(ch_vcf)
      if ( params.qc && ! params.popStrat) {
        sample_qc(buildConversion.out.bed, buildConversion.out.bim, ch_fam)
        variant_qc(sample_qc.out.bed, sample_qc.out.bim, sample_qc.out.fam)
      } else if ( params.qc && params.popStrat) {
        sample_qc(buildConversion.out.bed, buildConversion.out.bim, ch_fam)
        popStrat(sample_qc.out.bed, sample_qc.out.bim, sample_qc.out.fam)
        variant_qc(popStrat.out.bed, popStrat.out.bim, popStrat.out.fam)
      }
      if ( params.impute ) {
        imputation(variant_qc.out.bed, variant_qc.out.bim, variant_qc.out.fam)
	      postImputation(imputation.out.imputed, variant_qc.out.fam)
	      if (params.gwas) {
	    gwas(postImputation.out.bed, postImputation.out.bim, postImputation.out.fam, popStrat.out.covar)
	}
      } else if ( !params.impute && params.gwas ) {
         if (params.popStrat) {
	   gwas(variant_qc.out.bed, variant_qc.out.bim, variant_qc.out.fam, popStrat.out.covar)
	 } else if (!params.popStrat) {
     // use dummy covar
	   gwas(variant_qc.out.bed, variant_qc.out.bim, variant_qc.out.fam, dummy_covar)
        }
      }
    }

    // workflow without build conversion
    if ( !params.convertBuild ) {
      if ( params.qc && !params.popStrat ) {
        sample_qc(ch_bed, ch_bim, ch_fam)	
        variant_qc(sample_qc.out.bed, sample_qc.out.bim, sample_qc.out.fam)
      } else if ( params.qc && params.popStrat ) {
        sample_qc(ch_bed, ch_bim, ch_fam)
        popStrat(sample_qc.out.bed, sample_qc.out.bim, sample_qc.out.fam)
        variant_qc(popStrat.out.bed, popStrat.out.bim, popStrat.out.fam)  
      }
      // imputation 
      if ( params.impute ) {
        imputation(variant_qc.out.bed, variant_qc.out.bim, variant_qc.out.fam)
	      postImputation(imputation.out.imputed, variant_qc.out.fam)
	      if (params.gwas && params.popStrat) {
          gwas(postImputation.out.bed, postImputation.out.bim, postImputation.out.fam, popStrat.out.covar)
        } else if (params.gwas && ! params.popStrat) {
          gwas(postImputation.out.bed, postImputation.out.bim, postImputation.out.fam, dummy_covar)
        }
      } 
      if ( !params.impute && params.gwas ) {
        if (params.popStrat) {
          gwas(variant_qc.out.bed, variant_qc.out.bim, variant_qc.out.fam, popStrat.out.covar)
        } else if (!params.popStrat) {
          // use dummy covar
          gwas(variant_qc.out.bed, variant_qc.out.bim, variant_qc.out.fam, dummy_covar)
      }
    }
  }
}

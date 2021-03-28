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
include {preImputation} from './workflows/preImputation.nf'
include {gwas} from './workflows/gwas.nf'
include {download_core} from './workflows/download_db.nf'
include {download_impute} from './workflows/download_db.nf'

if (params.help) {
  printHelp()
  System.exit(0)
}

if (!params.convert_build && !params.qc && !params.popStrat && !params.impute && !params.gwas && !params.download_db && !params.pre_impute && !params.post_impute) {
  println("${params}")
  println("Please specify some workflow options")
  println("------------------------------------")
  printHelp()
  System.exit(1)
}

println """
=================================================================
snpQT is ready to make your single-nucleotide polymorphisms cute!
v1.0, February 2021
Parameters in effect:
${params}
=================================================================
        """.stripIndent()

// throw errors on invalid workflow combinations --------------------------
if (!params.download_db ==~ "core" || !params.download_db ==~ "impute") {
  println("Please use --download_db core or --download_db impute")
  System.exit(1)
}

if (params.convert_build) {
    if (!params.vcf) {
      println("Please supply a vcf.gz file for build conversion with --vcf")
      println("Use --help to print help")
      System.exit(1)
      }
  else if (!params.fam) {
    println("Please supply a .fam file for build conversion with --fam")
    println("Use --help to print help")
    System.exit(1)
    }
  } else if (!params.convert_build) {
    if (params.vcf) {
      println("--vcf only compatible with --convert_build")
      println("Please supply plink input files with --bed --bim --fam")
      println("Use --help to print help")
      System.exit(1)
      }
  }
  
if (params.qc) {
	 if (params.bed && !params.bim || !params.bed && params.bim ) {
		println("--bed and --bim must be supplied together")
		println("Use --help to print help")
		System.exit(1)
	  } else if (!params.bed || !params.bim || !params.bed ) {
		println("Missing --fam, --bed and --bim parameters")
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
	} else if (!params.qc && params.pre_impute) {
	  println("--pre_impute requires --qc")
	  println("Use --help to print help")
	  System.exit(1)
	} else if (!params.qc && params.gwas) {
	  println("--gwas requires --qc")
	  println("Please rerun with --qc")
	  println("Use --help to print help")
	  System.exit(1)
	}

if (params.impute && params.pre_impute) {
    println("--pre_impute is not combined with --impute")
    println("--impute supports Pre-Imputation and Post-Imputation automatically")
    println("Please rerun without --pre_impute")
	println("Use --help to print help")
    System.exit(1)
  }
  
if (params.gwas && params.pre_impute) {
    println("--pre_impute is not combined with --gwas")
    println("--pre_impute is designed to prepare your VCF for imputation in an external server")
    println("If you wish to run local imputation use --impute")
    println("Please rerun --gwas without --pre_impute")
	println("Use --help to print help")
    System.exit(1)
  }


if (params.convert_build || params.qc ) {
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
  if ( params.convert_build ) {
    Channel
      .fromPath(params.vcf, checkIfExists: true)
      .set{ ch_vcf }
	Channel
      .fromPath(params.fam, checkIfExists: true)
      .set{ ch_fam }
  }
  
  if ( !params.convert_build && params.qc ) {
    Channel
      .fromPath(params.bed, checkIfExists: true)
      .set{ ch_bed }
    Channel
      .fromPath(params.bim, checkIfExists: true)
      .set{ ch_bim }
    Channel
      .fromPath(params.fam, checkIfExists: true)
      .set{ ch_fam }
  } else if (params.convert_build && params.qc) {
      Channel
        .fromPath(params.fam, checkIfExists: true)
        .set{ ch_fam }
  }

  main:
    if ( params.download_db == "core" ) {
      download_core()
    } else if (params.download_db == "impute") {
      download_impute()
    } 

    // workflow with build conversion
    if ( params.convert_build) {
      buildConversion(ch_vcf, ch_fam)
      if ( params.qc && ! params.popStrat) {
        sample_qc(buildConversion.out.bed, buildConversion.out.bim, buildConversion.out.fam)
        variant_qc(sample_qc.out.bed, sample_qc.out.bim, sample_qc.out.fam)
      } else if ( params.qc && params.popStrat) {
        sample_qc(buildConversion.out.bed, buildConversion.out.bim, buildConversion.out.fam)
        popStrat(sample_qc.out.bed, sample_qc.out.bim, sample_qc.out.fam)
        variant_qc(popStrat.out.bed, popStrat.out.bim, popStrat.out.fam)
      }
	  // pre-imputation
	  if ( params.pre_impute ) {
        preImputation(variant_qc.out.bed, variant_qc.out.bim, variant_qc.out.fam)
		println("Go ahead and upload your VCF to an external Imputation server.")
      }
	  // imputation & GWAS
      if ( params.impute ) {
	    preImputation(variant_qc.out.bed, variant_qc.out.bim, variant_qc.out.fam)
        imputation(preImputation.out.vcf, preImputation.out.idx)
		postImputation(imputation.out.imputed, variant_qc.out.fam)
		if ( params.gwas ) {
			gwas(postImputation.out.bed, postImputation.out.bim, postImputation.out.fam, variant_qc.out.covar)
		}
      } else if ( !params.impute && params.gwas ) {
        gwas(variant_qc.out.bed, variant_qc.out.bim, variant_qc.out.fam, variant_qc.out.covar)
      }
    }

    // workflow without build conversion
    if ( !params.convert_build ) {
      if ( params.qc && !params.popStrat ) {
        sample_qc(ch_bed, ch_bim, ch_fam)	
        variant_qc(sample_qc.out.bed, sample_qc.out.bim, sample_qc.out.fam)
      } else if ( params.qc && params.popStrat ) {
        sample_qc(ch_bed, ch_bim, ch_fam)
        popStrat(sample_qc.out.bed, sample_qc.out.bim, sample_qc.out.fam)
        variant_qc(popStrat.out.bed, popStrat.out.bim, popStrat.out.fam)  
      }
	  // pre-imputation
	  if ( params.pre_impute ) {
        preImputation(variant_qc.out.bed, variant_qc.out.bim, variant_qc.out.fam)
		println("Go ahead and upload your VCF to an external Imputation server.")
      }
      // imputation 
      if ( params.impute ) {
        imputation(variant_qc.out.bed, variant_qc.out.bim, variant_qc.out.fam)
        postImputation(imputation.out.imputed, variant_qc.out.fam)
        if (params.gwas) {
          gwas(postImputation.out.bed, postImputation.out.bim, postImputation.out.fam, variant_qc.out.covar)
        } 
      } else if ( !params.impute && params.gwas ) {
        gwas(variant_qc.out.bed, variant_qc.out.bim, variant_qc.out.fam, variant_qc.out.covar)
      }
    }
}

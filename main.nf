#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {printHelp} from './modules/help.nf'

// import subworkflows
include {buildConversion} from './workflows/buildConversion.nf'
include {sample_qc} from './workflows/sample_qc.nf'
include {variant_qc} from './workflows/variant_qc.nf'
include {pop_strat} from './workflows/popStrat.nf'
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

if (!params.convert_build && !params.qc && !params.pop_strat && !params.impute && !params.gwas && !params.download_db && !params.pre_impute && !params.post_impute) {
  println """
  =================================================================
  snpQT is ready to make your single-nucleotide polymorphisms cute!
  v0.1.4 - Fluffy penguin, 2021-07-14
  Parameters in effect:
  ${params}
  =================================================================
	  """.stripIndent()
  println("Please specify some workflow options")
  println("------------------------------------")
  printHelp()
  System.exit(1)
}

println """
=================================================================
snpQT is ready to make your single-nucleotide polymorphisms cute!
v0.1.3 - Fluffy penguin, 2021-06-15
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
	} else if (!params.fam) {
		println("Please supply a .fam file for build conversion with --fam")
		println("Use --help to print help")
		System.exit(1)
	} else if (params.input_build != 37 && params.input_build != 38) {
		println("Build conversion workflow supports only build 37 and 38")
		println("Please use --input_build 37 or --input_build 38 and --output_build 38 or --output_build 37 for build conversion")
		println("Use --help to print help")
		System.exit(1)
	} else if (params.output_build != 37 && params.output_build != 38) {
		println("Build conversion workflow supports only build 37 and 38")
		println("Please use --input_build 37 and --output_build 38 or --input_build 38 and --output_build 37 for build conversion")
		println("Use --help to print help")
		System.exit(1)
	}
  } else if (!params.convert_build && !params.post_impute) {
    if (params.vcf) {
      println("--vcf only compatible with --convert_build")
      println("Please supply plink input files with --bed --bim --fam")
      println("Use --help to print help")
      System.exit(1)
      }
  }
  
if (params.qc) {
		if (!params.convert_build){
		 if (params.bed && !params.bim || !params.bed && params.bim ) {
			println("--bed and --bim must be supplied together")
			println("Use --help to print help")
			System.exit(1)
		  } else if (!params.bed || !params.bim || !params.bed ) {
			println("Missing --fam, --bed and --bim parameters")
			println("Use --help to print help")
			System.exit(1)
		  }
		}
	} else if (!params.qc && params.pop_strat) {
	  println("--pop_strat requires --qc")
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
	} else if (!params.qc && params.gwas && !params.post_impute) {
	  println("--gwas requires --qc or --post_impute")
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

//Pre-Imputation, Imputation and Post-Imputation compatibility errors

if ( params.impute && params.pre_impute ) {
    println("--pre_impute is not combined with --impute")
    println("--impute supports Pre-Imputation and Post-Imputation automatically")
    println("Please rerun keeping only --pre_impute or --impute workflow parameters")
	println("Use --help to print help")
    System.exit(1)
  }
  
if ( params.impute && params.post_impute ) {
    println("--post_impute is not combined with --impute")
    println("--impute supports Pre-Imputation and Post-Imputation automatically")
    println("Please rerun keeping only --post_impute or --impute workflow parameters")
	println("Use --help to print help")
    System.exit(1)
  }
  
if ( params.pre_impute && params.post_impute ) {
    println("--post_impute is not combined with --pre_impute")
    println("Please rerun keeping only --post_impute or --pre_impute workflow parameters")
	println("Use --help to print help")
    System.exit(1)
  }
  
if ( params.gwas && params.pre_impute ) {
    println("--pre_impute is not combined with --gwas")
    println("--pre_impute is designed to prepare your VCF for imputation in an external server")
    println("If you wish to run local imputation use --impute")
    println("Please rerun --gwas without --pre_impute")
	println("Use --help to print help")
    System.exit(1)
  }

if ( params.post_impute && params.qc){
	println("--post_impute is not combined with --qc")
    println("--post_impute expects a VCF file and a .fam file whereas --qc expects binary plink files")
    println("Please rerun --post_impute without --qc")
	println("Use --help to print help")
    System.exit(1)	
}

if ( params.post_impute && params.gwas){
	println("--post_impute is not combined with --gwas")
    println("If you wish to run --qc please first run --post-impute alone, and then use the output binary plink files to run other snpQT workflows like --qc, --pop_strat and --gwas")
    println("Please rerun --post_impute providing only a --vcf and a --fam file ")
	println("Use --help to print help")
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
  
  if (params.post_impute) {
	Channel
	  .fromPath(params.fam, checkIfExists: true)
	  .set{ ch_fam }
	Channel
	  .fromPath(params.vcf, checkIfExists: true)
	  .set{ ch_imp }
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
      if ( params.qc && ! params.pop_strat) {
        sample_qc(buildConversion.out.bed, buildConversion.out.bim, buildConversion.out.fam)
        variant_qc(sample_qc.out.bed, sample_qc.out.bim, sample_qc.out.fam)
      } else if ( params.qc && params.pop_strat) {
        sample_qc(buildConversion.out.bed, buildConversion.out.bim, buildConversion.out.fam)
        pop_strat(sample_qc.out.bed, sample_qc.out.bim, sample_qc.out.fam)
        variant_qc(pop_strat.out.bed, pop_strat.out.bim, pop_strat.out.fam)
      }
	  // pre-imputation
	  if ( params.pre_impute ) {
        preImputation(variant_qc.out.bed, variant_qc.out.bim, variant_qc.out.fam)
		println("Go ahead and upload your VCF to an external Imputation server.")
      }
	  // imputation & GWAS
      if ( params.impute ) {
	    preImputation(variant_qc.out.bed, variant_qc.out.bim, variant_qc.out.fam)
            imputation(preImputation.out.vcf)
		postImputation(imputation.out.imputed_vcf, variant_qc.out.fam)
		if ( params.gwas ) {
			gwas(postImputation.out.bed, postImputation.out.bim, postImputation.out.fam, variant_qc.out.covar)
		}
      } else if ( !params.impute && params.gwas ) {
        gwas(variant_qc.out.bed, variant_qc.out.bim, variant_qc.out.fam, variant_qc.out.covar)
      }
    }

    // workflow without build conversion
    if ( !params.convert_build ) {
      if ( params.qc && !params.pop_strat ) {
        sample_qc(ch_bed, ch_bim, ch_fam)	
        variant_qc(sample_qc.out.bed, sample_qc.out.bim, sample_qc.out.fam)
      } else if ( params.qc && params.pop_strat ) {
        sample_qc(ch_bed, ch_bim, ch_fam)
        pop_strat(sample_qc.out.bed, sample_qc.out.bim, sample_qc.out.fam)
        variant_qc(pop_strat.out.bed, pop_strat.out.bim, pop_strat.out.fam)  
      }
	  // pre-imputation without imputation
	  if ( params.pre_impute ) {
        preImputation(variant_qc.out.bed, variant_qc.out.bim, variant_qc.out.fam)
      }
      // local imputation 
      if ( params.impute ) {
	    preImputation(variant_qc.out.bed, variant_qc.out.bim, variant_qc.out.fam)
        imputation(preImputation.out.vcf)
        postImputation(imputation.out.imputed_vcf, variant_qc.out.fam)
        if ( params.gwas ) {
          gwas(postImputation.out.bed, postImputation.out.bim, postImputation.out.fam, variant_qc.out.covar)
        } 
      } else if ( !params.impute && params.gwas && params.qc ) {
        gwas(variant_qc.out.bed, variant_qc.out.bim, variant_qc.out.fam, variant_qc.out.covar)
      }
    }
	
	// post-imputation workflow
	if ( params.post_impute ) {
        postImputation(ch_imp, ch_fam) 
	}
}

params.inbed = "../../results/sample_qc/sample_variant_qc.*"
params.inbim = "../../results/sample_qc/sample_variant_qc.bim"
params.infam = "../../results/sample_qc/sample_variant_qc.fam"
params.outdir = "$baseDir/../../results"
outdir = params.outdir + '/popStrat/'

log.info """\
         snpQT step03: Imputation
         input file: ${params.inbim}
         outdir: ${params.outdir}
         """
         .stripIndent()

Channel.fromPath( params.inbed ).into { inbed; inbed_snpflip } 
Channel.fromPath( params.inbim ).into { inbim; inbim_snpflip } 
Channel.fromPath( params.infam ).into { infam; infam_snpflip } 

Channel.fromPath("$SNPQT_DB_DIR/human_g1k_v37.fasta").into{ g37_D1 ; g37_D8 }
Channel.fromPath("$SNPQT_DB_DIR/All_20180423.vcf.gz").set{ dbSNP }
Channel.fromPath("$SNPQT_DB_DIR/All_20180423.vcf.gz.tbi").set{ dbSNP_index }

// Pre-imputation 
// =============================================================================

// STEP D1: Check for strand issues, ambiguous and duplicated SNPs using snpflip
process run_snpflip {
    container 'snpflip'

    input:
    file g37_D1
    file inbed_snpflip
    file inbim_snpflip
    file infam_snpflip

    output:
    file "D1*" into D1

    shell:
    '''
    snpflip -b sample_variant_qc.bim \
        -f !{g37_D1} \
        -o D1
    '''

}

// STEP D2: Remove ambiguous SNPs ---------------------------------------------
// STEP D3: Remove one of each pair of duplicated SNPs 
// STEP D4: Flip all SNPs that are on the reverse strand 

process flip_snps {
    input:
    file inbed
    file inbim
    file infam 
    file D1

    output:
    file 'D4*' into D4

    shell:
    '''
    # D2: ambiguous 
    plink -bfile sample_variant_qc \
        --exclude D1.ambiguous \
        --make-bed \
        --out D2 

    # D3: duplicates
    plink --bfile D2 \
        --list-duplicate-vars ids-only suppress-first
    
    plink --bfile D2 \
        --exclude plink.dupvar \
        --make-bed \
        --out D3

    # D4: reverse 
    plink --bfile D3 \
        --flip D1.reverse \
        --make-bed \
        --out D4  
    '''
}

// STEP D5: Convert Plink file into VCF ---------------------------------------
// STEP D6: bgzip VCF and then VCF needs to be indexed/sorted 
// STEP D7: Convert .vcf.gz file to .bcf file ---------------------------------

process to_bcf {
    input:
    file D4 

    output:
    file 'D7.bcf' into D7

    shell:
    '''
    plink --bfile D4 \
        --recode vcf bgz \
        --keep-allele-order \
        --out D6
    bcftools convert -Ou D6.vcf.gz > D7.bcf
    '''
}

// STEP D8: Check and fix the REF allele --------------------------------------

process check_ref_allele {
    input:
    file D7 
    file dbSNP
    file dbSNP_index
    file g37_D8

    output:
    file 'D8.bcf' into D8

    shell:
    '''
    bcftools +fixref !{D7} \
        -Ob -o D8.bcf -- \
        -d -f !{g37_D8} \
        -i !{dbSNP}
    '''
}

// STEP D9: Sort the BCF ------------------------------------------------------
// STEP D10: Convert .bcf file to .vcf.gz file --------------------------------
// STEP D11: Index the vcf.gz -------------------------------------------------

process sort_to_vcf {
    input:
    file D8

    output:
    file 'D11.vcf.gz' into D11
    file 'D11.vcf.gz.csi' into D11_index
    
    shell:
    '''
    bcftools sort !{D8} | bcftools convert -Oz > D11.vcf.gz
    bcftools index D11.vcf.gz 
    '''
}

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
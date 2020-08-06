params.inbed = "../../results/sample_qc/sample_variant_qc.*"
params.inbim = "../../results/sample_qc/sample_variant_qc.bim"
params.infam = "../../results/sample_qc/sample_variant_qc.fam"
params.outdir = "$baseDir/../../results"
outdir = params.outdir + '/imputation/'

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
Channel.fromPath("$SNPQT_DB_DIR/genetic_maps.b37.tar.gz").set{ shapeit4_map }
Channel
    .fromFilePairs( ["$SNPQT_DB_DIR/thousand_genomes/*.vcf.gz", "$SNPQT_DB_DIR/thousand_genomes/*.vcf.gz.tbi"] )
    .set{ thousand_genomes } 

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
Channel.from(1..22).set{ split_list } // groovy range 

process split_chrom {
    input:
    file D11
    file D11_index
    each chr from split_list 

    output:
    tuple val(chr), file('D12.vcf.gz') into D12
    tuple val(chr), file('D12.vcf.gz.csi') into D12_index 

    shell:
    '''
    bcftools view -r !{chr} !{D11} -Oz -o D12.vcf.gz
    bcftools index D12.vcf.gz
    '''
}

// STEP D14: Perform phasing using shapeit4 -----------------------------------
// join vcf file channel with index file channel based on chrom value 
// then combine so each tuple element has a shapeit4 map file 
// kind of like 'each', but for tuples
D12_combined = D12.join(D12_index).combine(shapeit4_map)

process phasing {
    container 'shapeit4' 

    input:
    tuple chr, 'D12.vcf.gz', 'D12.vcf.gz.csi', 'genetic_maps.b37.tar.gz' from D12_combined 

    output:
    tuple chr, 'D14.vcf.gz' into D14, D14_idx

    shell:
    '''
    tar -xzf genetic_maps.b37.tar.gz
    gunzip chr!{chr}.b37.gmap.gz # decompress the chromosome we need 
    
    shapeit4 --input D12.vcf.gz \
        --map chr!{chr}.b37.gmap \
        --region !{chr} \
        --thread 1 \
        --output D14.vcf.gz \
        --log log_chr.txt     
    '''
}

// STEP D15: Index phased chromosomes -----------------------------------------

process index {
    input:
    tuple chr, 'D14.vcf.gz' from D14_idx

    output:
    tuple chr, 'D14.vcf.gz.csi' into D14_csi
    
    shell:
    '''
    bcftools index D14.vcf.gz
    '''
}

// Imputation 
// =============================================================================

// Note: STEP D16 is taken care of by Dockerfile 
// STEP D17: Convert vcf reference genome into a .imp5 format for each chromosome

process imp5convert {
    container 'impute5'

    input:
    tuple val(x), path(chrom) from thousand_genomes

    output:
    tuple val(CHR), '1k_b37_reference_chr.imp5' into D17  

    shell:
    '''
    # CHR is chromosome integer, x is filename prefix from file pair 
    CHR=`echo !{x} | cut -d '.' -f 2 | tr -dc '0-9'`
    # chrom is a file pair of vcf and idx, idx not useful after imp5
    CHROM=`echo !{chrom} | cut -d ' ' -f 1` # there must be a nicer way >:(     

    imp5Converter --h $CHROM \
        --r $CHR \
        --o 1k_b37_reference_chr.imp5
    '''
}

// STEP D18: Perform imputation using impute5 ---------------------------------

// Finished!
log.info """\
         snpQT: make your SNPs cute 
         install_impute: Process reference data downloaded by download_impute.sh
         input directory : ${params.indir}
         """
         .stripIndent()

Channel.fromPath(params.indir + "/*genotypes.vcf.gz").set{ ref_vcf } 
Channel.fromPath(params.indir + "/human_g1k_v37.fasta").set{ h37 }
// Step F3

process qc {
    publishDir params.indir, mode: 'copy', overwrite: true
    
    input:
    each vcf from ref_vcf
    file h37

    output:
    file "*updated.vcf.gz"
    file "*updated.vcf.gz.*"
    
    shell:
    '''
    bcftools norm -m-any --check-ref w -f !{h37} !{vcf} | \
		bcftools norm -Oz --rm-dup both -o !{vcf.baseName}_updated.vcf.gz
    bcftools index !{vcf.baseName}_updated.vcf.gz
    '''
}
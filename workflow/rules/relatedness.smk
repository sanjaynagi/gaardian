rule ngsRelate:
    """
    Run NGSRelate on VCF files
    """
    input:
        vcf = "resources/vcfs/{dataset}_{chrom}.biallelic.vcf.gz"
    output:
        "results/relatedness/ngsRelate_{chrom}"
    log:
        log = "logs/ngsRelate/{chrom}.log"
    params:
        tag = 'GT',
        basedir=workflow.basedir,
    shell:
        """
        {params.basedir}/scripts/NgsRelate/ngsRelate -h {input.vcf} -O {output} -c 1 -T {params.tag} 2> {log}
        """
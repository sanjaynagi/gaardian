
rule ZarrToVCF:
    input:
        ZarrGeno = "resources/snp_genotypes/all/{_set}/{chrom}/calldata/GT/",
        ZarrSiteFilters = "resources/site_filters/dt_20200416/{sitefilter}/{chrom}/variants/filter_pass",
        ZarrPOS = "resources/snp_genotypes/all/sites/{chrom}/variants/POS/",
        metadata = "resources/metadata/samples.meta.csv"
    output:
        multiallelicVCF = "resources/vcfs/{dataset}_{chrom}.vcf",
        biallelicVCF = "resources/vcfs/{dataset}_{chrom}.biallelic.vcf""
    params:
        basedir=workflow.basedir
    script:
        "{params.basedir}/scripts/ZarrToVCF.py"


rule bgzip:
    input:
        calls = "resources/vcfs/{dataset}_{chrom}{allelism}.vcf",
    output:
        calls_gz = "resources/vcfs/{dataset}_{chrom}{allelism}.vcf.gz",
    log:
        "logs/bgzip/{chrom}_{allelism}.log",
    shell:
        """
        bgzip {input.calls} 2> {log}
        """

rule BCFtoolsIndex:
    input:
        calls = "resources/vcfs/{dataset}_{chrom}{allelism}.vcf.gz",
    output:
        calls_gz = "resources/vcfs/{dataset}_{chrom}{allelism}.vcf.gz.csi",
    log:
        "logs/bcftoolsIndex/{chrom}_{allelism}.log",
    shell:
        """
        bcftools index {input.calls} 2> {log}
        """

rule tabix:
    input:
        calls = "resources/vcfs/{dataset}_{chrom}{allelism}.vcf.gz",
    output:
        calls_tbi = "resources/vcfs/{dataset}_{chrom}{allelism}.vcf.gz.tbi",
    log:
        "logs/tabix/{chrom}_{allelism}.log",
    shell:
        """
        tabix {input.calls} 2> {log}
        """


    
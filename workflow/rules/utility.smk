

rule GenomeIndex:
    """
    Index the reference genome with samtools
    """
    input:
        ref=config["reference"]["genome"],
    output:
        idx=config["reference"]["genome"] + ".fai",
    log:
        "logs/GenomeIndex.log",
    wrapper:
        "v0.69.0/bio/samtools/faidx"


rule ZarrToVCF:
    """
    Write out biallelic and multiallelic VCF files from provided Zarr files 
    """
    input:
        genotypes = config['Zarr']['Genotypes'],
        siteFilters = config['Zarr']['SiteFilters'],
        Positions = config['Zarr']['Positions'],
    output:
        multiallelicVCF = "resources/vcfs/{dataset}_{chrom}.multiallelic.vcf",
        biallelicVCF = "resources/vcfs/{dataset}_{chrom}.biallelic.vcf"
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        basedir=workflow.basedir,
        metadata = config['metadata'],
    script:
        "{params.basedir}/scripts/ZarrToVCF.py"



gzippedVCF = getVCFs(gz=True, bothAllelisms=True)
rule BGZip:
    """
    This is overwriting log files at the
    """
    input:
        calls = getVCFs(gz=False, bothAllelisms=True)
    output:
        calls_gz = gzippedVCF
    log:
        "logs/bgzip/{chrom}_{allelism}.log" if config['VCF']['activate'] is False else "logs/bgzip/{chrom}.log"
    shell:
        """
        bgzip {input.calls} 2> {log}
        """

rule BcftoolsIndex:
    input:
        calls = getVCFs(gz=True)
    output:
        calls_gz = "resources/vcfs/{dataset}_{chrom}.{allelism}.vcf.gz.csi",
    log:
        "logs/bcftoolsIndex/{dataset}_{chrom}.{allelism}.log",
    shell:
        """
        bcftools index {input.calls} 2> {log}
        """

rule Tabix:
    input:
        calls = getVCFs(gz=True)
    output:
        calls_tbi = "resources/vcfs/{dataset}_{chrom}.{allelism}.vcf.gz.tbi",
    log:
        "logs/tabix/{dataset}_{chrom}_{allelism}.log",
    shell:
        """
        tabix {input.calls} 2> {log}
        """
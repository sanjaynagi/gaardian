

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
        genotypes = getZarrArray(type_="Genotypes") if not cloud else [],
        positions = getZarrArray(type_='Positions') if not cloud else [],
        siteFilters = getZarrArray(type_ = "SiteFilters") if not cloud else [],
        refPath = config['Zarr']['REFPath'] if not cloud else [],
        altPath = config['Zarr']['ALTPath'] if not cloud else []
    output:
        multiallelicVCF = expand("resources/vcfs/{dataset}_{{chrom}}.multiallelic.vcf", dataset=dataset),
        biallelicVCF = expand("resources/vcfs/{dataset}_{{chrom}}.biallelic.vcf", dataset=dataset)
    conda:
        "../envs/pythonGenomics.yaml"
    log:
        "logs/ZarrToVCF/{chrom}.log"
    params:
        cloud=cloud,
        ag3_sample_sets = ag3_sample_sets,
        basedir=workflow.basedir,
        metadata = config['metadata'],
        dataset = dataset
    script:
        "../scripts/ZarrToVCF.py"




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
        "logs/bgzip/{contig}_{allelism}.log" if config['VCF']['activate'] is False else "logs/bgzip/{contig}.log"
    shell:
        """
        bgzip {input.calls} 2> {log}
        """

rule BcftoolsIndex:
    input:
        calls = getVCFs(gz=True)
    output:
        calls_gz = "resources/vcfs/{dataset}_{contig}.{allelism}.vcf.gz.csi",
    log:
        "logs/bcftoolsIndex/{dataset}_{contig}.{allelism}.log",
    shell:
        """
        bcftools index {input.calls} 2> {log}
        """

rule Tabix:
    input:
        calls = getVCFs(gz=True)
    output:
        calls_tbi = "resources/vcfs/{dataset}_{contig}.{allelism}.vcf.gz.tbi",
    log:
        "logs/tabix/{dataset}_{contig}_{allelism}.log",
    shell:
        """
        tabix {input.calls} 2> {log}
        """
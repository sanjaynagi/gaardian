

rule pca:
    """
    Perform principal components analysis 
    """
    input:
        genotypes = getZarrArray(type_="Genotypes"),
        positions = getZarrArray(type_='Positions'),
        siteFilters = getZarrArray(type_ = "SiteFilters"),
    output:
        htmlAll = expand("results/PCA/{dataset}.{{contig}}.html", dataset = dataset),
        pngAll = expand("results/PCA/{dataset}.{{contig}}.png", dataset = dataset),
        html = expand("results/PCA/{cohort}.{{contig}}.html", cohort=PCAcohorts['cohortNoSpaceText']),
        png = expand("results/PCA/{cohort}.{{contig}}.png", cohort=PCAcohorts['cohortNoSpaceText']),
    log:
        log = "logs/pca/{contig}.log"
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        metadata = config['metadata'],
        dataset = config['dataset'],
        data = "results/PCA/data",
    script:
        "../scripts/pca.py"



rule f2HapLength:
    """
    Find lengths of haplotypes
    """
    input:
        f2variantPairs = "results/f2variantPairs.tsv"
    output:
        "results/f2HapLengths_{contig}.tsv"
    log:
        log = "logs/f2variants/hapLength_{contig}.log"
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        metadata = config['metadata'],
        genotypes = getZarrArray(type_="Genotypes"),
        positions = getZarrArray(type_='Positions'),
    script:
        "../scripts/f2HaplotypeLength.py"




rule ngsRelate:
    """
    Run NGSRelate on VCF files
    """
    input:
        vcf = getVCFs(allelism='biallelic')
    output:
        "results/relatedness/ngsRelate.{dataset}.{contig}"
    log:
        log = "logs/ngsRelate/{dataset}_{contig}.log"
    params:
        tag = 'GT',
        basedir=workflow.basedir,
    threads: 24
    shell:
        """
        ../scripts/NgsRelate/ngsRelate -h {input.vcf} -O {output} -c 1 -T {params.tag} -p {threads} 2> {log}
        """
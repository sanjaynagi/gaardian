

rule pca:
    """
    Perform principal components analysis 
    """
    input:
        genotypes = getZarrArray(type_="Genotypes"),
        positions = getZarrArray(type_='Positions'),
        siteFilters = getZarrArray(type_ = "SiteFilters"),
    output:
        htmlAll = expand("results/PCA/{dataset}.{{chrom}}.html", dataset = dataset),
        pngAll = expand("results/PCA/{dataset}.{{chrom}}.png", dataset = dataset),
        html = expand("results/PCA/{cohort}.{{chrom}}.html", cohort=PCAcohorts['cohortNoSpaceText']),
        png = expand("results/PCA/{cohort}.{{chrom}}.png", cohort=PCAcohorts['cohortNoSpaceText']),
    log:
        log = "logs/pca/{chrom}.log"
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        metadata = config['metadata'],
        dataset = config['dataset'],
        data = "results/PCA/data",
    script:
        "../scripts/pca.py"


rule ngsRelate:
    """
    Run NGSRelate on VCF files
    """
    input:
        vcf = getVCFs(allelism='biallelic')
    output:
        "results/relatedness/ngsRelate.{dataset}.{chrom}"
    log:
        log = "logs/ngsRelate/{dataset}_{chrom}.log"
    params:
        tag = 'GT',
        basedir=workflow.basedir,
    threads: 24
    shell:
        """
        ../scripts/NgsRelate/ngsRelate -h {input.vcf} -O {output} -c 1 -T {params.tag} -p {threads} 2> {log}
        """
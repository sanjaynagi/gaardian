

rule pca:
    """
    Perform principal components analysis 
    """
    input:
        genotypes = getZarrArray(type_="Genotypes", cloud=cloud),
        positions = getZarrArray(type_='Positions', cloud=cloud),
        siteFilters = getZarrArray(type_ = "SiteFilters", cloud=cloud),
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
        cohortColumn = config['PopulationStructure']['PCA']['colourColumns'],
        cloud = cloud,
        ag3_sample_sets = ag3_sample_sets
    script:
        "../scripts/pca.py"


rule f2VariantLocate:
    """
    Find lengths of haplotypes
    """
    input:
        genotypes = getZarrArray(type_="Genotypes", cloud=cloud),
        positions = getZarrArray(type_='Positions', cloud=cloud),
        sitefilters = getZarrArray(type_='SiteFilters', cloud=cloud),
    output:
        f2variantPairs = "results/f2variantPairs.tsv"
    log:
        log = "logs/f2variants/f2VariantLocate.log"
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        metadata = config['metadata'],
        genotypes = getZarrArray(type_="Genotypes", cloud=cloud),
        positions = getZarrArray(type_='Positions', cloud=cloud),
        cloud = cloud,
        contigs = contigs,
        ag3_sample_sets = ag3_sample_sets
    script:
        "../scripts/f2VariantLocate.py"



rule f2HapLength:
    """
    Find lengths of haplotypes
    """
    input:
        genotypes = getZarrArray(type_="Genotypes", cloud=cloud),
        positions = getZarrArray(type_='Positions', cloud=cloud),
        f2variantPairs = "results/f2variantPairs.tsv"
    output:
        "results/f2HapLengths_{contig}.tsv"
    log:
        log = "logs/f2variants/hapLength_{contig}.log"
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        metadata = config['metadata'],
        genotypes = getZarrArray(type_="Genotypes", cloud=cloud),
        positions = getZarrArray(type_='Positions', cloud=cloud),
        cloud = cloud,
        ag3_sample_sets = ag3_sample_sets
    script:
        "../scripts/f2HaplotypeLength.py"




rule ngsRelate:
    """
    Run NGSRelate on VCF files
    """
    input:
        vcf = "resources/vcfs/wholegenome/{dataset}.biallelic.vcf.gz",
        csi = "resources/vcfs/wholegenome/{dataset}.biallelic.vcf.gz.csi",
        tbi = "resources/vcfs/wholegenome/{dataset}.biallelic.vcf.gz.tbi",
    output:
        "results/relatedness/ngsRelate.{dataset}"
    log:
        log = "logs/ngsRelate/{dataset}.log"
    params:
        tag = 'GT',
        basedir=workflow.basedir,
    threads: 24
    shell:
        """
        {params.basedir}/scripts/NgsRelate/ngsRelate -h {input.vcf} -O {output} -c 1 -T {params.tag} -p {threads} 2> {log}
        """
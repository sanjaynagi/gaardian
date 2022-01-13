
## Optionally specifify PBS
## G123/

rule VariantsOfInterest:
    """
    This rule reports and plots allele frequencies of Variants of Interest specified in VariantsOfInterest.tsv.
    """
    input:
 #       genotypes = 
    output:
        "results/VariantsOfInterest/{}"
    log:
        "logs/VariantsOfInterest.log"
    params:
    script:
        "../scripts/VariantsOfInterest.py"


rule G12:
    """
    This rule performs G12 selection scans on each specified population
    """
    input:
        genotypes = getZarrArray(type_="Genotypes"),
        positions = getZarrArray(type_='Positions'),
        siteFilters = getZarrArray(type_ = "SiteFilters")
    output:
        plot = expand("results/selection/G12/G12_{cohort}.{{chrom}}.png", cohort=cohorts['cohortNoSpaceText']),
        tsv = expand("results/selection/G12/G12_{cohort}.{{chrom}}.tsv", cohort=cohorts['cohortNoSpaceText'])
    log:
        "logs/selection/G12.{chrom}.log"
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        metadata = config['metadata'],
        columns = config['metadataCohortColumns'],
        GarudsStat = "G12",
        windowSize = config['Selection']['G12']['windowSize'],
        windowStep = config['Selection']['G12']['windowStep'],
        cutHeight = config['Selection']['G12']['cutHeight'],
        minPopSize = 15
    script:
        "../scripts/GarudsStatistics.py"

rule G123:
    """
    This rule performs G123 selection scans on each specified population
    """
    input:
        genotypes = getZarrArray(type_="Genotypes"),
        positions = getZarrArray(type_='Positions'),
        siteFilters = getZarrArray(type_ = "SiteFilters")
    output:
        plot = expand("results/selection/G123/G123_{cohort}.{{chrom}}.png", cohort=cohorts['cohortNoSpaceText']),
        tsv = expand("results/selection/G123/G123_{cohort}.{{chrom}}.tsv", cohort=cohorts['cohortNoSpaceText'])
    log:
        "logs/selection/G123.{chrom}.log"
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        metadata = config['metadata'],
        columns = config['metadataCohortColumns'],
        GarudsStat = "G123",
        windowSize = config['Selection']['G123']['windowSize'],
        windowStep = config['Selection']['G123']['windowStep'],
        cutHeight = config['Selection']['G123']['cutHeight'],
        minPopSize = 15
    script:
        "../scripts/GarudsStatistics.py"


rule H12:
    """
    This rule performs H12 selection scans on each specified population
    """
    input:
        haplotypes = getZarrArray(type_="Genotypes"),
        positions = getZarrArray(type_='Positions'),
        siteFilters = getZarrArray(type_ = "SiteFilters")
    output:
        plot = expand("results/selection/H12/H12_{cohort}.{{chrom}}.png", cohort=cohorts['cohortNoSpaceText']),
        tsv = expand("results/selection/H12/H12_{cohort}.{{chrom}}.tsv", cohort=cohorts['cohortNoSpaceText'])
    log:
        "logs/selection/H12.{chrom}.log"
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        metadata = config['metadata'],
        columns = config['metadataCohortColumns'],
        windowSize = config['Selection']['H12']['windowSize'],
        windowStep = config['Selection']['H12']['windowStep'],
        minPopSize = 15
    script:
        "../scripts/GarudsStatistics.py"


rule PopulationBranchStatistic:
    """
    This rule performs PBS selection scans on each specified population
    """
    input:
        genotypes = getZarrArray(type_="Genotypes"),
        positions = getZarrArray(type_='Positions'),
        siteFilters = getZarrArray(type_ = "SiteFilters")
    output:
        plot = expand("results/selection/PBS/PBS_{cohort}.{{chrom}}.png", cohort=cohorts['cohortNoSpaceText']),
        tsv = expand("results/selection/PBS/PBS_{cohort}.{{chrom}}.tsv", cohort=cohorts['cohortNoSpaceText'])
    log:
        "logs/selection/PBS.{chrom}.log"
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        metadata = config['metadata'],
        columns = config['metadataCohortColumns'],
        comparatorColumn = config['Selection']['PBS']['metadataComparatorColumn'],
        windowSize = config['Selection']['PBS']['windowSize'],
        windowStep = config['Selection']['PBS']['windowStep'],
        minPopSize = 15
    script:
        "../scripts/PopulationBranchStatistic.py"



#from tools import get_colour_dict
#from tools import getCohorts
#import matplotlib.pyplot as plt
#import matplotlib
#import numpy as np


def getCohorts(metadata, columns=['species_gambiae_coluzzii', 'location'], comparatorColumn=None, minPopSize=15):
    
    # subset metadata dataFrame and find combinations of variables with more than minPopSize individuals
    cohorts = metadata[columns]
    cohorts = cohorts.groupby(columns).size().reset_index().rename(columns={0:'size'})
    cohorts = cohorts[cohorts['size'] > minPopSize][columns]
    
    if comparatorColumn != None:
        cols = [i for i in columns if i != comparatorColumn]
    else:
        cols = columns

    idxs = []
    for _, row in cohorts.iterrows():   
        # create the pandas metadata query for each cohort
        mycohortQuery = " & ".join([col + " == " + "'" + row.astype(str)[col] + "'" for col in cohorts.columns])
        # get indices of individuals in each cohort
        idxs.append(metadata.query(mycohortQuery).index.tolist())
    
    cohorts['indices'] = idxs
    cohorts['cohortText'] = cohorts[cols].agg(' | '.join, axis=1)
    cohorts['cohortNoSpaceText'] = cohorts['cohortText'].str.replace("|", ".", regex=False).str.replace(" ", "",regex=False)
    #colours = get_colour_dict(cohorts['species_gambiae_coluzzii'], palette="Set1")
    #cohorts['colour'] = cohorts['species_gambiae_coluzzii'].map(colours)
    if comparatorColumn != None: 
  #      cols = columns + ['cohortText', 'cohortNoSpaceText']
        cohorts = cohorts.pivot(index='cohortNoSpaceText', columns=comparatorColumn)
        return(cohorts.reset_index())

    return(cohorts.reset_index(drop=True))



def getZarrArray(type_="Genotype", all_contigs=False):

    if config['Zarr']['activate'] == True:
        Array = config['Zarr'][type_]
        if all_contigs == True:
            Array = Array.replace("{chrom}", "{{chrom}}")
            print(Array)
    elif config['Zarr']['activate'] == False:
        if type_ == "Genotype":
            Array = "resources/Zarr/{dataset}/{chrom}/calldata/GT" 
        elif type_ == "Positions":
            Array = "resources/Zarr/{dataset}/{chrom}/variants/POS"
        elif type_ == 'SiteFilters' and siteFilters is not None:
            Array = "resources/Zarr/{dataset}/{chrom}/variants/siteFilter"

        if all_contigs == True:
            Array = Array.replace("{chrom}", "{{chrom}}")
    else:
        Array = []

#    print(Array)
    return(Array)


def getVCFs(gz=True,allelism = 'biallelic', bothAllelisms=False):

    if config['VCF']['activate'] == True:
        genotypes = config['VCF'][allelism]
    elif gz == True:
        genotypes = expand("resources/vcfs/{dataset}_{{chrom}}.{{allelism}}.vcf.gz", dataset=config['dataset']) if bothAllelisms == True else expand("resources/vcfs/{dataset}_{{chrom}}.{allelism}.vcf.gz", dataset=config['dataset'], allelism=allelism)
    elif gz == False:
        genotypes = expand("resources/vcfs/{dataset}_{{chrom}}.{{allelism}}.vcf", dataset=config['dataset']) if bothAllelisms == True else expand("resources/vcfs/{dataset}_{{chrom}}.{allelism}.vcf", dataset=config['dataset'], allelism=allelism)

    return(genotypes)

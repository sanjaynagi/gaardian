#!/usr/bin/env python
# coding: utf-8

"""
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

import probetools as probe
import numpy as np
import pandas as pd
import allel
import dask.array as da
import scipy
import seaborn as sns
import matplotlib.pyplot as plt


# Garuds Selection Scans # 
cloud = snakemake.params['cloud']
ag3_sample_sets = snakemake.params['ag3_sample_sets']
contig = snakemake.wildcards['contig']
stat = snakemake.params['GarudsStat']
windowSize = snakemake.params['windowSize']
windowStep = snakemake.params['windowStep']

if not cloud:
    haplotypePath = snakemake.input['haplotypes'] if stat in ['H1', 'H12', 'H2/1'] else []
    positionsPath = snakemake.input['positions']
else:
    haplotypePath = []
    positionsPath = []

# Load metadata 
if cloud:
    import malariagen_data
    ag3 = malariagen_data.Ag3(pre=True)
    metadata = ag3.sample_metadata(sample_sets=ag3_sample_sets)
else:
    metadata = pd.read_csv(snakemake.params['metadata'], sep="\t")

# Load arrays 
haps, pos = probe.loadZarrArrays(haplotypePath, 
        positionsPath, 
        siteFilterPath=None, 
        haplotypes=True, 
        cloud=cloud, 
        contig=contig)


# Define functions
def clusterMultiLocusGenotypes(gnalt, cut_height=0.1, metric='euclidean', g=2):
    """
    Hierarchically clusters genotypes and calculates G12 statistic. 
    """
    # cluster the genotypes in the window
    dist = scipy.spatial.distance.pdist(gnalt.T, metric=metric)
    if metric in {'hamming', 'jaccard'}:
        # convert distance to number of SNPs, easier to interpret
        dist *= gnalt.shape[0]

    Z = scipy.cluster.hierarchy.linkage(dist, method='single')
    cut = scipy.cluster.hierarchy.cut_tree(Z, height=cut_height)[:, 0]
    cluster_sizes = np.bincount(cut)
    #clusters = [np.nonzero(cut == i)[0] for i in range(cut.max() + 1)] #returns indices of individuals in each cluster
    
    # get freq of clusters and sort by largest freq
    freqs = cluster_sizes/gnalt.shape[1]
    freqs = np.sort(freqs)[::-1]
    
    # calculate garuds statistic
    gStat = np.sum(freqs[:g])**2 + np.sum(freqs[g:]**2)
    
    return(gStat)


def garudsStat(stat, geno, pos, cut_height=None, metric='euclidean', window_size=1200, step_size=600):
    
    """
    Calculates G12/G123/H12
    """
        
    # Do we want to cluster the Multi-locus genotypes (MLGs), or just group MLGs if they are identical
    if stat == "G12":
        garudsStat = allel.moving_statistic(geno, clusterMultiLocusGenotypes, size=window_size, step=step_size, metric=metric, cut_height=cut_height, g=2)
    elif stat == "G123":
        garudsStat = allel.moving_statistic(geno, clusterMultiLocusGenotypes, size=window_size, step=step_size, metric=metric, cut_height=cut_height, g=3)
    elif stat == "H12":
        garudsStat,_,_,_ = allel.moving_garud_h(geno, size=window_size, step=step_size)
    else:
        raise ValueError("Statistic is not G12/G123/H12")

    midpoint = allel.moving_statistic(pos, np.median, size=window_size, step=step_size)
    
    return(garudsStat, midpoint)



#### Load cohort data and their indices in genotype data
### run garudStat for that query. already loaded contigs 

cohorts = probe.getCohorts(metadata=metadata, 
                    columns=snakemake.params.columns, 
                    minPopSize=snakemake.params.minPopSize, excludepath="resources/sib_group_table.csv")


# Loop through each cohort, manipulate genotype arrays and calculate chosen Garuds Statistic
for idx, cohort in cohorts.iterrows():

    # get indices for haplotype Array and filter
    hapInds = np.sort(np.concatenate([np.array(cohort['indices'])*2, np.array(cohort['indices']*2)+1]))
    gt_cohort = haps.take(hapInds, axis=1)

    probe.log(f"--------- Running {stat} on {cohort['cohortText']} | Chromosome {contig} ----------")
    probe.log("filter to biallelic segregating sites")

    ac_cohort = gt_cohort.count_alleles(max_allele=3)
    # N.B., if going to use to_n_alt later, need to make sure sites are 
    # biallelic and one of the alleles is the reference allele
    loc_sites = ac_cohort.is_biallelic()
    gt_seg = gt_cohort.compress(loc_sites, axis=0)
    pos_seg = pos[loc_sites]

    # probe.log(f"compute input data for {stat}")
    pos_seg = pos_seg

    # calculate G12/G123/H12 and plot figs 
    gStat, midpoint = garudsStat(stat=stat,
                                geno=gt_seg, 
                                pos=pos_seg, 
                                cut_height=None,
                                metric='euclidean',
                                window_size=windowSize,
                                step_size=windowStep)

    probe.windowedPlot(statName=stat, 
                cohortText = cohort['cohortText'],
                cohortNoSpaceText= cohort['cohortNoSpaceText'],
                values=gStat, 
                midpoints=midpoint,
                prefix=f"results/selection/{stat}", 
                contig=contig,
                colour=cohort['colour'],
                ymin=0,
                ymax=0.5)

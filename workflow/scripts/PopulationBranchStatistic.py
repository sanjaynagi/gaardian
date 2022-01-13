#!/usr/bin/env python
# coding: utf-8

"""
PBS

TODO

scripts still needs - way of getting comparator column order and ensuring only 2 options (ie. case control or dead alive)
                    - needs outgroup path input
                    - need to make GAARD metadata and run
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

from tools import loadZarrArrays, getCohorts, windowedPlot, log
from pathlib import Path
import numpy as np
import pandas as pd
import allel
import dask.array as da
import seaborn as sns
import matplotlib.pyplot as plt


# PBS Selection Scans # 
chrom = snakemake.wildcards['chrom']
stat = "PBS"
windowSize = snakemake.params['windowSize']
windowStep = snakemake.params['windowStep']
genotypePath = snakemake.input['genotypes']
positionsPath = snakemake.input['positions']
siteFilterPath = snakemake.input['siteFilters']

# Outgroup data
outgroupPath = f"/home/sanj/ag1000g/data/phase3/snp_genotypes/all/AG1000G-ML-B/{chrom}/calldata/GT/"
outgroupMetaPath = "/home/sanj/ag1000g/data/phase3/metadata/general/AG1000G-ML-B/samples.meta.csv"
Mali2004Meta = pd.read_csv(outgroupMetaPath)
species = pd.read_csv("/home/sanj/ag1000g/data/phase3/metadata/species_calls_20200422/AG1000G-ML-B/samples.species_aim.csv")
Mali2004Meta = Mali2004Meta.merge(species)

# Read metadata 
metadata = pd.read_csv(snakemake.params['metadata'], sep=",")
metadata['location'] = metadata['location'].str.split(".").str.get(0)


# Load arrays
snps, pos = loadZarrArrays(genotypePath, positionsPath, siteFilterPath=siteFilterPath, haplotypes=False)

### Load outgroup Arrays and subset to each species, storing
snpsOutgroup, pos = loadZarrArrays(outgroupPath, positionsPath, siteFilterPath=siteFilterPath, haplotypes=False)
snpsOutgroupDict = {}

for sp in ['gambiae', 'coluzzii']:
    sp_bool = Mali2004Meta['species_gambiae_coluzzii'] == sp
    snpsOutgroupDict[sp] =  snpsOutgroup.compress(sp_bool, axis=1)


#### Load cohort data and their indices in genotype data
### run garudStat for that query. already loaded chroms 

cohorts = getCohorts(metadata=metadata, 
                    columns=snakemake.params.columns, 
                    comparatorColumn=snakemake.params.comparatorColumn,
                    minPopSize=snakemake.params.minPopSize)





# Loop through each cohort, manipulate genotype arrays and calculate chosen Garuds Statistic
for idx, cohort in cohorts.iterrows():

    log(f"--------- Running {stat} on {cohort['cohortText'].to_list()} | Chromosome {chrom} ----------")
    log("filter to biallelic segregating sites")    
    species = cohort['species_gambiae_coluzzii'].to_list()[0]

    if len(cohort['indices']['alive']) < snakemake.params.minPopSize:
        continue
    elif len(cohort['indices']['dead']) < snakemake.params.minPopSize:
        continue

    ac_cohort = snps.count_alleles(max_allele=3).compute()
    # N.B., if going to use to_n_alt later, need to make sure sites are 
    # biallelic and one of the alleles is the reference allele
    ref_ac = ac_cohort[:, 0]
    loc_sites = ac_cohort.is_biallelic() & (ref_ac > 0)
    gt_seg = da.compress(loc_sites, snps, axis=0)
    pos_seg = da.compress(loc_sites, pos, axis=0)
    
    log(f"compute input data for {stat}")
    pos_seg = pos_seg.compute()
    
    ac_out = allel.GenotypeArray(da.compress(loc_sites, snpsOutgroupDict[species], axis=0)).count_alleles()
    ac_alive = allel.GenotypeArray(gt_seg).take(cohort['indices']['alive'], axis=1).count_alleles()
    ac_dead = allel.GenotypeArray(gt_seg).take(cohort['indices']['dead'], axis=1).count_alleles()
        
    assert ac_out.shape[0] == pos_seg.shape[0], "Array Outgroup/POS are the wrong length"
    assert ac_alive.shape[0] == pos_seg.shape[0], "Array alive/POS are the wrong length"
    assert ac_dead.shape[0] == pos_seg.shape[0], "Arrays dead/POS the wrong length"

    log("calculate PBS and plot figs")
    # calculate PBS and plot figs 
    pbsArray = allel.pbs(ac_alive, ac_dead, ac_out, 
                window_size=size, window_step=step, normed=True)
    midpoint = allel.moving_statistic(pos_seg, np.mean, size=size, step=step)
    
    windowedPlot(statName=stat, 
                cohortText = cohort['cohortText'].to_numpy()[0],
                cohortNoSpaceText= cohort['cohortNoSpaceText'].to_numpy()[0],
                values=pbsArray, 
                midpoints=midpoint,
                prefix=f"results/selection/{stat}", 
                chrom=chrom,
                colour=cohort['colour'].to_numpy()[0],
                ymin=-0.5,
                ymax=0.5,
                save=True)

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
contig = snakemake.wildcards['contig']
stat = "PBS"
windowSize = snakemake.params['windowSize']
windowStep = snakemake.params['windowStep']
genotypePath = snakemake.input['genotypes']
positionsPath = snakemake.input['positions']
siteFilterPath = snakemake.input['siteFilters']

# Outgroup data
outgroupPath = snakemake.input['outgroupPath']
outgroupMetaPath = snakemake.input['outgroupMetaPath']
Mali2004Meta = pd.read_csv(outgroupMetaPath)
species = pd.read_csv("/home/sanj/ag1000g/data/phase3/metadata/species_calls_20200422/AG1000G-ML-B/samples.species_aim.csv")
Mali2004Meta = Mali2004Meta.merge(species)

# Read metadata 
metadata = pd.read_csv(snakemake.params['metadata'], sep="\t")

# Load arrays
snps, pos = loadZarrArrays(genotypePath, positionsPath, siteFilterPath=siteFilterPath)

### Load outgroup Arrays and subset to each species, storing
snpsOutgroup, pos = loadZarrArrays(outgroupPath, positionsPath, siteFilterPath=siteFilterPath)
snpsOutgroupDict = {}

for sp in ['gambiae', 'coluzzii']:
    sp_bool = Mali2004Meta['species_gambiae_coluzzii'] == sp
    snpsOutgroupDict[sp] =  snpsOutgroup.compress(sp_bool, axis=1)

#### Load cohort data and their indices in genotype data
### run garudStat for that query. already loaded contigs
cohorts = getCohorts(metadata=metadata,
                    columns=snakemake.params.columns,
                    comparatorColumn=snakemake.params.comparatorColumn,
                    minPopSize=snakemake.params.minPopSize)
cohorts = cohorts.dropna()

# Get name for phenotype of interest
pheno1, pheno2 = cohorts['indices'].columns.to_list()

# Loop through each cohort, manipulate genotype arrays and calculate chosen Garuds Statistic
for idx, cohort in cohorts.iterrows():

    log(f"--------- Running {stat} on {cohort['cohortText'].to_list()} | Chromosome {contig} ----------")
    log("filter to biallelic segregating sites")    
    species = cohort['species_gambiae_coluzzii'].to_list()[0]
    if len(cohort['indices'][pheno1]) < snakemake.params.minPopSize:
        continue
    elif len(cohort['indices'][pheno2]) < snakemake.params.minPopSize:
        continue
    elif cohort['indices'][pheno1] == 'NaN':
        continue
    elif cohort['indices'][pheno2] == 'NaN':
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
    ac_pheno1 = allel.GenotypeArray(gt_seg).take(cohort['indices'][pheno1], axis=1).count_alleles()
    ac_pheno2 = allel.GenotypeArray(gt_seg).take(cohort['indices'][pheno2], axis=1).count_alleles()
        
    assert ac_out.shape[0] == pos_seg.shape[0], "Array Outgroup/POS are the wrong length"
    assert ac_pheno1.shape[0] == pos_seg.shape[0], "Array phenotype1/POS are the wrong length"
    assert ac_pheno2.shape[0] == pos_seg.shape[0], "Arrays phenotype2/POS the wrong length"

    log("calculate PBS and plot figs")
    # calculate PBS and plot figs 
    pbsArray = allel.pbs(ac_pheno1, ac_pheno2, ac_out, 
                window_size=windowSize, window_step=windowStep, normed=True)
    midpoint = allel.moving_statistic(pos_seg, np.mean, size=windowSize, step=windowStep)
    
    windowedPlot(statName=stat, 
                cohortText = cohort['cohortText'].to_numpy()[0],
                cohortNoSpaceText= cohort['cohortNoSpaceText'].to_numpy()[0],
                values=pbsArray, 
                midpoints=midpoint,
                prefix=f"results/selection/{stat}", 
                contig=contig,
                colour=cohort['colour'].to_numpy()[0],
                ymin=-0.3,
                ymax=0.3,
                save=True)

#!/usr/bin/env python
# coding: utf-8

"""
Zarr to VCF. Needs providing with REF and ALT and allpositions files.
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

from tools import loadZarrArrays, log
from pathlib import Path
import numpy as np
import zarr
import pandas as pd
import allel
import dask.array as da
from datetime import date
from pathlib import Path


# Zarr to VCF # 
cloud = snakemake.params['cloud']
ag3_sample_sets = snakemake.params['ag3_sample_sets'] if cloud else []
contig = snakemake.wildcards['contig']
dataset = snakemake.params['dataset']
genotypePath = snakemake.input['genotypes'] if not cloud else []
positionsPath = snakemake.input['positions'] if not cloud else []
siteFilterPath = snakemake.input['siteFilters'] if not cloud else []
refPath = snakemake.input['refPath']
altPath = snakemake.input['altPath']

# Load metadata 
if cloud:
    import malariagen_data
    ag3 = malariagen_data.Ag3()
    metadata = ag3.sample_metadata(sample_sets=ag3_sample_sets)
else:
    metadata = pd.read_csv(snakemake.params['metadata'], sep="\t")

def write_vcf_header(vcf_file, contig):
    """
    Writes a VCF header.
    """
    
    print('##fileformat=VCFv4.1', file=vcf_file)
    # write today's date
    today = date.today().strftime('%Y%m%d')
    print('##fileDate=%s' % today, file=vcf_file)
    # write source
    print('##source=scikit-allel-%s + ZarrToVCF.py' % allel.__version__, file=vcf_file)
    #write refs and contigs 
    print('##reference=resources/reference/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa', file=vcf_file)
    print('##contig=<ID=2R,length=61545105>', file=vcf_file) if contig == '2R' else None
    print('##contig=<ID=3R,length=53200684>', file=vcf_file) if contig == '3R' else None 
    print('##contig=<ID=2L,length=49364325>', file=vcf_file) if contig == '2L' else None
    print('##contig=<ID=3L,length=41963435>', file=vcf_file) if contig == '3L' else None
    print('##contig=<ID=X,length=24393108>', file=vcf_file) if contig == 'X' else None
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=vcf_file)


def ZarrToPandasToVCF(vcf_file, genotypePath, positionsPath, siteFilterPath, contig, nchunks=50, snpfilter = "segregating"):
    
    """
    Converts genotype and POS arrays to vcf, using pd dataframes in chunks. 
    Segregating sites only. Needs REF and ALT arrays.
    """
    
    #if file exists ignore and skip
    myfile = Path(f"{vcf_file}.gz")
    if myfile.is_file():
        print(f"File {vcf_file}.gz Exists...")
        return
    
    log(f"Loading array for {contig}...")

    geno, pos = loadZarrArrays(genotypePath, positionsPath, siteFilterPath=siteFilterPath, cloud=cloud, contig=contig, sample_sets=ag3_sample_sets, haplotypes=False)
    allpos = allel.SortedIndex(zarr.open_array(positionsPath)[:])
    ref_alt_filter = allpos.locate_intersection(pos)[0]
    
    refs = zarr.open_array(refPath.format(contig=contig))[:][ref_alt_filter]
    alts = zarr.open_array(altPath.format(contig=contig))[:][ref_alt_filter]
    
    if snpfilter == "segregating":
        log("Find segregating sites...")
        flt = geno.count_alleles().is_segregating()
        geno = geno.compress(flt, axis=0)
        positions = pos[flt]
        refs = refs[flt].astype(str)
        alts = [a +"," + b + "," + c for a,b,c in alts[flt].astype(str)]
    elif snpfilter == 'biallelic01':
        log("Finding biallelic01 sites...")
        flt = geno.count_alleles().is_biallelic_01()
        geno = geno.compress(flt, axis=0)
        positions = pos[flt]
        refs = refs[flt].astype(str)
        alts = [a for a,b,c in alts[flt].astype(str)]
    else:
        assert np.isin(snpfilter, ['segregating', "biallelic01"]).any(), "incorrect snpfilter value"
  
    log("calculating chunks sizes...")
    chunks = np.round(np.arange(0, geno.shape[0], geno.shape[0]/nchunks)).astype(int).tolist()
    chunks.append(geno.shape[0])

    for idx, chunk in enumerate(chunks[:-1]):

        gn = geno[chunks[idx]:chunks[idx+1]].compute()
        pos = positions[chunks[idx]:chunks[idx+1]]
        ref = refs[chunks[idx]:chunks[idx+1]]
        alt = alts[chunks[idx]:chunks[idx+1]]
        
        # Contruct SNP info DF
        vcf_df = pd.DataFrame({'#CHROM': contig,
                 'POS': pos,
                 'ID': '.',
                 'REF': ref,
                 'ALT': alt,
                 'QUAL': '.',
                 'FILTER': '.',
                 'INFO':'.',
                'FORMAT': 'GT'})

        log(f"Pandas SNP info DataFrame constructed...{idx}")

        # Geno to VCF
        vcf = pd.DataFrame(gn.to_gt().astype(str), columns=metadata['specimen'])
        log("Concatenating info and genotype dataframes...")
        vcf = pd.concat([vcf_df, vcf], axis=1)

        log(f"Pandas Genotype data constructed...{idx}")

        if (idx==0) is True:
            with open(f"{vcf_file}", 'w') as vcfheader:
                    write_vcf_header(vcfheader, contig)

        log("Writing to .vcf")

        vcf.to_csv(vcf_file, 
                   sep="\t", 
                   index=False,
                   mode='a',
                  header=(idx==0), 
                  line_terminator="\n")



### MAIN ####

#contigs = ['2L', '2R', '3L', '3R', 'X']
#ZarrToPandasToVCF(f"../resources/vcfs/ag3_gaardian_{contig}.multiallelic.vcf", genotypePath, positionsPath, siteFilterPath, contig, snpfilter="segregating")


ZarrToPandasToVCF(f"resources/vcfs/{dataset}_{contig}.biallelic.vcf", genotypePath, positionsPath, siteFilterPath, contig, snpfilter="biallelic01")

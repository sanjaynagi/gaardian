#!/usr/bin/env python
# coding: utf-8


"""
Zarr to VCF. Needs providing with REF and ALT and allpositions files.
"""

import sys
#sys.stderr = open(snakemake.log[0], "w")

import numpy as np
import zarr
import pandas as pd
import allel
import dask.array as da
from datetime import date
from pathlib import Path

def ZarrToPandasToHaplotypeVCF(vcf_file, metadata, sample_sets, contig, sample_query=None, analysis='gamb_colu', nchunks=50, sampleNameColumn = 'partner_sample_id'):
    
    """
    Converts genotype and POS arrays to vcf, using pd dataframes in chunks. 
    Segregating sites only. Needs REF and ALT arrays.
    """
    
    #if file exists ignore and skip
    myfile = Path(f"{vcf_file}.gz")
    if myfile.is_file():
        print(f"File {vcf_file}.gz Exists...")
        return
    
    print(f"Loading array for {contig}...")

    ds_haps = ag3.haplotypes(contig, sample_sets=sample_sets, sample_query=sample_query, analysis=analysis)
    sample_ids = ds_haps['sample_id'].values
    metadata = metadata.set_index('sample_id').loc[sample_ids, :].reset_index()
    positions = ds_haps['variant_position']
    geno = allel.GenotypeDaskArray(ds_haps['call_genotype'])

    refs = ds_haps['variant_allele'][:,0].compute().values.astype(str)
    alts = ds_haps['variant_allele'][:,1].compute().values.astype(str)
  
    print("calculating chunks sizes...")
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

        print(f"Pandas SNP info DataFrame constructed...{idx}")

        # Geno to VCF
        vcf = pd.DataFrame(np.char.replace(gn.to_gt().astype(str), "/", "|"), columns=metadata[sampleNameColumn])
        print("Concatenating info and genotype dataframes...")
        vcf = pd.concat([vcf_df, vcf], axis=1)

        print(f"Pandas Genotype data constructed...{idx}")

        if (idx==0) is True:
            with open(f"{vcf_file}", 'w') as vcfheader:
                    write_vcf_header(vcfheader, contig)

        print("Writing to .vcf")

        vcf.to_csv(vcf_file, 
                   sep="\t", 
                   index=False,
                   mode='a',
                  header=(idx==0), 
                  lineterminator="\n")

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

import sys


# Zarr to VCF # 
ag3_sample_sets = ["1288-VO-UG-DONNELLY-VMF00168","1288-VO-UG-DONNELLY-VMF00219"] #'1244-VO-GH-YAWSON-VMF00149'  #snakemake.params['ag3_sample_sets'] if cloud else []
contig = sys.argv[1] #'2L' #snakemake.wildcards['contig']
dataset = 'llineup' #snakemake.params['dataset']
sampleNameColumn = 'partner_sample_id'
sample_query = 'taxon == "gambiae"'

import malariagen_data
ag3 = malariagen_data.Ag3(pre=True)
metadata = ag3.sample_metadata(sample_sets=ag3_sample_sets, sample_query=sample_query)

print(f"Running for {contig}...")

### MAIN ####
ZarrToPandasToHaplotypeVCF(
     f"resources/vcfs/{dataset}_{contig}.haplotypes.vcf", 
     metadata=metadata,
     sample_query=sample_query,
     contig=contig, 
     nchunks=20, 
     sample_sets=ag3_sample_sets
    )
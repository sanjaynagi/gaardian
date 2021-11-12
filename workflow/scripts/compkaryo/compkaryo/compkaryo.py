#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 15 16:57:04 2018

@author: beccalove
"""

import argparse
import pkgutil
import sys

import numpy as np
import allel


## read in the predictive SNPs for the inversion of interest

## extract the predictive SNPs from the supplied callset

## calculate the average genotype across the predictive SNPs

## return the average genotype and the # of sites across which it was calculated

## please note: you can modify this software to work in other systems

## or with other SNPs by changing inversionDict


inversionDict = {"2La" : ("2L", "targets/2La_targets.txt"),
                 "2Rj" : ("2R", "targets/2Rj_targets.txt"),
                 "2Rb" : ("2R", "targets/2Rb_targets.txt"),
                 "2Rc_col" : ("2R", "targets/2Rc_col_targets.txt"),
                 "2Rc_gam" : ("2R", "targets/2Rc_gam_targets.txt"),
                 "2Rd" : ("2R", "targets/2Rd_targets.txt"),
                 "2Ru" : ("2R", "targets/2Ru_targets.txt"),
                 "test" : ("fooY", "targets/test_targets_clean.txt")}

def parse_args(custom_args=None):

    '''Interpret command-line arguments.'''

    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="path to variant call format file \
                containing the genotypes")
    parser.add_argument("inversion", help="inversion to be classified",
                        choices=["2La", "2Rj", "2Rb", "2Rc_col", "2Rc_gam",
                                 "2Rd", "2Ru", "test"])
    parser.add_argument("-o", "--out",
                        help="name of the results file")
    parser.add_argument("-s", "--samples",
                        help=\
                        "samples to include; file or space-separated list",
                        nargs='+')
    parser.add_argument("-t", "--totals",
                        help="total # sites supporting each genotype",
                        action='store_true')
    parser.add_argument("-p", "--ploidy", help="Ploidy level in VCF", default=2)

    args = parser.parse_args(custom_args)

    return args

def get_numbers_dict(ploidy):

    """ 
    Initialise dict for sci-kit allel.read_vcf() 
    Allows for variable ploidy in vcf (useful for pool/RNA-Seq)
    """
    numbers = {
        'samples':1,
        'variants/CHROM': 1,
        'variants/POS': 1,
        'variants/ID': 1,
        'variants/REF': 1,
        'variants/ALT': 'A',
        'variants/QUAL': 1,
        'variants/DP': 1,
        'variants/AN': 1,
        'variants/AC': 'A',
        'variants/AF': 'A',
        'variants/MQ': 1,
        'variants/ANN': 1,
        'calldata/DP': 1,
        'calldata/GT': ploidy,
        'calldata/GQ': 1,
        'calldata/HQ': 2,
        'calldata/AD': 'R',
        'calldata/MQ0': 1,
        'calldata/MQ': 1,
        }
        
    return numbers

def import_data(callset_path):

    '''Read in the VCF in the appropriate format.'''

    numbers = get_numbers_dict(args.ploidy)

    callset = allel.read_vcf(callset_path,
                             numbers=numbers,
                             fields=['samples', 'calldata/GT',
                                     'variants/CHROM', 'variants/FILTER',
                                     'variants/POS', 'variants/REF',
                                     'variants/ALT'])

    return callset

def import_inversion(inversion):

    '''Load the tag SNPs for the appropriate inversion.'''

    path = inversionDict[inversion][1]

    targets_raw = pkgutil.get_data("compkaryo", path)

    targets = np.array([int(entry)\
                        for entry in targets_raw.decode().splitlines()\
                        if not entry.startswith('#')])

    return targets

def extract_vtbl_indices(targets, callset, chrom=None):

    '''Find where in the genotypes the tag SNPs are located.'''

    if not chrom:
        chrom = inversionDict[args.inversion][0]

    indices = []

    for site in targets:

        where = np.where( (callset["variants/POS"] == site) &\
                         (callset["variants/CHROM"] == chrom))

        if len(where[0]) > 0:

            indices.append(where[0][0])

    return indices

def create_samples_bool(callset):

    '''Take only a subset of samples.'''

    if len(args.samples) == 1 and args.samples[0].endswith('.txt'):

        samples_file_handle = args.samples[0]

        samples = np.genfromtxt(samples_file_handle, dtype="str")

    else:

        samples = np.array(args.samples)

    samples_bool =\
    np.array([sample in samples for sample in callset["samples"]])

    return samples_bool

def calculate_genotype_at_concordant_sites(callset, indices, samples_bool=None,
                                           totals=False):

    '''Calculate the average number of alternate alleles for each specimen at
    each tag SNP.'''

    genos = allel.GenotypeArray(callset["calldata/GT"]).subset(sel0=indices)

    if samples_bool is not None:

        genos = genos.subset(sel1=samples_bool)

    alt_count = genos.to_n_alt()

    is_called = genos.is_called()

    av_gts = np.mean(np.ma.MaskedArray(alt_count,
                                       mask=~is_called), axis=0).data

    if totals:

        match_dict = {0: None, 1: None, 2: None}

        for value in [0, 1, 2]:

            n_matches = np.sum(np.ma.MaskedArray(alt_count,
                                                 mask=~is_called) == value,
                                axis=0).data

            match_dict[value] = n_matches

    total_sites = np.sum(is_called, axis=0)

    if totals:

        data = av_gts, total_sites, match_dict[0], match_dict[1], match_dict[2]

    else:

        data = av_gts, total_sites

    return data

def main(args):

    '''Extract tag SNPs and desired specimens and calculate the average
    number of alternate alleles.'''

    callset = import_data(args.vcf)

    target_list = import_inversion(args.inversion)

    indices = extract_vtbl_indices(target_list, callset)

    samples_bool = None

    if args.samples:

        samples_bool = create_samples_bool(callset)

    #bi_bool = extract_biallelic_SNPs(callset, indices)

    if args.totals:

        av_gts, total_sites, num_0, num_1, num_2 =\
        calculate_genotype_at_concordant_sites(callset, indices, samples_bool,
                                               args.totals)

    else:

        av_gts, total_sites =\
        calculate_genotype_at_concordant_sites(callset, indices, samples_bool,
                                               args.totals)

    if not len(av_gts) == len(total_sites):

        raise ValueError("mean genotypes and total sites differ in length")

    if args.out is None:

        out = sys.stdout

    else:

        out = open(args.out, 'w')

    try:

        for i, gts in enumerate(av_gts):

            if args.totals:

                record = (str(gts), str(total_sites[i]), str(num_0[i]),
                          str(num_1[i]), str(num_2[i]))

            else:

                record = (str(gts), str(total_sites[i]))

            out_record = "\t".join(record) + "\n"

            out.write(out_record)

    finally:

        if args.out is not None:

            out.close()

if __name__ == "__main__":

    args = parse_args()
    main(args)
    
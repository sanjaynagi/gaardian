#!/usr/bin/env python
# coding: utf-8
"""
"""

import numpy as np
import pandas as pd
import malariagen_data
import allel
from malariagen_data.util import CacheMiss

def h12_gwss(
        contig,
        analysis,
        biallelic_mask,
        window_size,
        sample_sets=None,
        sample_query=None,
        cohort_size=30,
        random_seed=42,
    ):
        """Run h12 GWSS.

        Parameters
        ----------
        contig: str
            Chromosome arm (e.g., "2L")
        analysis : {"arab", "gamb_colu", "gamb_colu_arab"}
            Which phasing analysis to use. If analysing only An. arabiensis, the
            "arab" analysis is best. If analysing only An. gambiae and An.
            coluzzii, the "gamb_colu" analysis is best. Otherwise, use the
            "gamb_colu_arab" analysis.
        window_size : int
            The size of windows used to calculate h12 over.
        sample_sets : str or list of str, optional
            Can be a sample set identifier (e.g., "AG1000G-AO") or a list of
            sample set identifiers (e.g., ["AG1000G-BF-A", "AG1000G-BF-B"]) or a
            release identifier (e.g., "3.0") or a list of release identifiers.
        sample_query : str, optional
            A pandas query string which will be evaluated against the sample
            metadata e.g., "taxon == 'coluzzii' and country == 'Burkina Faso'".
        cohort_size : int, optional
            If provided, randomly down-sample to the given cohort size.
        random_seed : int, optional
            Random seed used for down-sampling.

        Returns
        -------
        x : numpy.ndarray
            An array containing the window centre point genomic positions.
        h12 : numpy.ndarray
            An array with h12 statistic values for each window.

        """
        # change this name if you ever change the behaviour of this function, to
        # invalidate any previously cached data
        name = "ag3_h12_gwss_v1"

        params = dict(
            contig=contig,
            analysis=analysis,
            biallelic_mask=biallelic_mask,
            window_size=window_size,
            sample_sets=ag3._prep_sample_sets_arg(sample_sets=sample_sets),
            sample_query=sample_query,
            cohort_size=cohort_size,
            random_seed=random_seed,
        )

        try:
            results = ag3.results_cache_get(name=name, params=params)

        except CacheMiss:
            results = _h12_gwss(**params)
            ag3.results_cache_set(name=name, params=params, results=results)

        x = results["x"]
        h12 = results["h12"]

        return x, h12

def _h12_gwss(
        contig,
        biallelic_mask,
        analysis,
        window_size,
        sample_sets,
        sample_query,
        cohort_size,
        random_seed,
    ):

        ds_haps = ag3.haplotypes(
            region=contig,
            analysis=analysis,
            sample_query=sample_query,
            sample_sets=sample_sets,
            cohort_size=cohort_size,
            random_seed=random_seed,
        )

        print(f"removing biallelics, {biallelic_mask.sum()} snps, array is shape: {ds_haps['call_genotype'].shape}")
        ds_haps = ds_haps.isel(variants=biallelic_mask)
        print(f"after, array is shape: {ds_haps['call_genotype'].shape}")
 
        gt = allel.GenotypeDaskArray(ds_haps["call_genotype"].data)

        with ag3._dask_progress(desc="Load haplotypes"):
            ht = gt.to_haplotypes().compute()
        pos = ds_haps["variant_position"].values

        h1, h12, h123, h2_h1 = allel.moving_garud_h(ht, size=window_size)

        x = allel.moving_statistic(pos, statistic=np.mean, size=window_size)

        results = dict(x=x, h12=h12)

        return results


def get_biallelic_mask(partner_sample_ids, cohort, contig, analysis):
    """
    Get a bialleic boolean mask which we will pass to H12_gwss() and apply to each phenotype/randomisation.
    Different boolean mask for each cohort.
    """

    print(f"Getting biallelic mask for {cohort}")
    ds_haps = ag3.haplotypes(
        region=contig,
        analysis=analysis,
        sample_query=f"partner_sample_id in {partner_sample_ids}",
        sample_sets=None,
        cohort_size=None,
        random_seed=None,
    )

    gt = allel.GenotypeDaskArray(ds_haps["call_genotype"].data)

    ac = gt.count_alleles()
    bial_bool = ac.is_biallelic()
    return(bial_bool)

# H12 Selection Scans # 
contig = snakemake.wildcards['contig']
cohort = snakemake.wildcards['cohort']
stat = "H12"

ag3 = malariagen_data.Ag3(pre=True)
metadata = ag3.sample_metadata("3.2")
meta = metadata.query("sex_call == 'F'")
sibs = pd.read_csv("resources/sib_group_table.csv", sep="\t")
exclude = sibs.query("keep == False")['sample.name']
metadata = meta.query("partner_sample_id not in @exclude").reset_index(drop=True)

randomisations = pd.read_csv("phenotype_randomisations.csv", sep="\t")
#### make sure all metadata samples are in the randomisations file (so GAARD data)
metadata = metadata.query(f"partner_sample_id in {randomisations['specimen'].to_list()}")
ids = metadata['partner_sample_id']
# make sure randomisations have only samples in the metadata, so excluding sibs and females
random = randomisations.query("specimen in @ids")


cohort_ids = random.query(f"population == @cohort")['specimen'].to_list()
biallelic_mask = get_biallelic_mask(cohort_ids, contig, cohort=cohort, analysis='gamb_colu')


alive_dead_dict = {}
alive_dead_dict['alive'] = random.query(f"phenotype == 'alive' & population == @cohort")['specimen']
alive_dead_dict['dead']  = random.query(f"phenotype == 'dead' & population == @cohort")['specimen']

for pheno, sample_names in alive_dead_dict.items():
# Loop through each cohort, manipulate genotype arrays and calculate chosen Garuds Statistic

    print(f"--------- Running {stat} on {cohort} {pheno} | Chromosome {contig} ----------")
    x, h12 = h12_gwss(
        contig=contig, 
        analysis='gamb_colu',
        biallelic_mask=biallelic_mask,
        window_size=1000,
        sample_query=f"partner_sample_id in {alive_dead_dict[pheno].to_list()}",
        cohort_size=None
    )

    h12_df = pd.DataFrame({'x':x, 'h12':h12})
    h12_df.to_csv(f"results/selection/H12/H12_{cohort}.{pheno}.{contig}.tsv", sep="\t")
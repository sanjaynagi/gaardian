### Probe
___


Probe is a snakemake workflow to apply a range of population genomic analyses to whole genome sequencing data, in any diploid organism. 
Data can be supplied in VCF or Zarr format, or *An. gambiae* s.l data can be accessed directly in the cloud through the [malariagen_data API](https://malariagen.github.io/vector-data/ag3/api.html). 
Called SNPs should be of high confidence, although a site filter mask can additionally be applied.

The workflow is a WIP, most modules are functional but nothing can be guaranteed :P 

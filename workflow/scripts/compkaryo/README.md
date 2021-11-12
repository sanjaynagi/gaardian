# How to use compkaryo

Compkaryo (short for "computational karyotyping") is a Python package for taking a VCF file generated from arm 2R or 2L of _Anopheles coluzzii_ and/or _gambiae_, and returning predictions for karyotypes of up to six inversions on those arms. Available inversions are: 2La; 2Rb; 2Rc; 2Ru; and for _An. gambiae_, 2Rd and 2Rj. The method behind these predictions is described in our publication, <a href=https://doi.org/10.1534/g3.119.400445>Love et al 2019</a> DOI: https://doi.org/10.1534/g3.119.400445.

Compkaryo does not currently do any filtering of the VCF. The file should be filtered to mask low-quality genotypes (we recommend masking everything with a GQ less than 20) before running compkaryo. Including low-quality genotypes in the calculation may __significantly__ skew the results.

Compkaryo requires Python 3 and depends on numpy and scikit-allel.

__Usage:__

python compkaryo.py vcf inversion

__Options:__

__-h, --help__      Print the help message

__-o, --out__      Specify an output file. If no output file is given, output prints to standard out.

__-s, --samples__   Specify sample names to include. Names can be given as a text file, with one name per line, or as a space-separated list on the command line.

__-t, --totals__    Output the total number of sites matching each possible inversion genotype (0, 1, 2).

 

Compkaryo currently operates on one inversion at a time. Use a loop to analyze multiple inversions.

GT16ErrorModel distribution
===========================
GT16ErrorModel(Double **epsilon**, Double **delta**, Alignment **alignment**)
-----------------------------------------------------------------------------

An error model distribution on an phased genotype alignment.

### Parameters

- Double **epsilon** - the sequencing and amplification error probability.
- Double **delta** - the allelic drop out probability.
- Alignment **alignment** - the genotype alignment.

### Return type

- Alignment

### Reference

Alexey Kozlov, Joao Alves, Alexandros Stamatakis, David Posada (2021). CellPhy: accurate and fast probabilistic inference of single-cell phylogenies from scDNA-seq data. bioRxiv 2020.07.31.230292.[https://doi.org/10.1101/2020.07.31.230292](https://doi.org/10.1101/2020.07.31.230292)


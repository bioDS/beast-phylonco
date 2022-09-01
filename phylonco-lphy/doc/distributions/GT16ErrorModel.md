GT16ErrorModel distribution
===========================
GT16ErrorModel([Double](../types/Double.md) **epsilon**, [Double](../types/Double.md) **delta**, [Alignment](../types/Alignment.md) **alignment**)
--------------------------------------------------------------------------------------------------------------------------------------------------

An error model distribution on an phased genotype alignment.

### Parameters

- [Double](../types/Double.md) **epsilon** - the sequencing and amplification error probability.
- [Double](../types/Double.md) **delta** - the allelic drop out probability.
- [Alignment](../types/Alignment.md) **alignment** - the genotype alignment.

### Return type

[Alignment](../types/Alignment.md)

### Reference

Alexey Kozlov, Joao Alves, Alexandros Stamatakis, David Posada (2021). CellPhy: accurate and fast probabilistic inference of single-cell phylogenies from scDNA-seq data. bioRxiv 2020.07.31.230292.[https://doi.org/10.1101/2020.07.31.230292](https://doi.org/10.1101/2020.07.31.230292)


# Phylonco
This is a BEAST 2 package for Bayesian inference of molecular data for cancer evolution. This package implements error models and substitution models for inference of timed trees in a Bayesian MCMC framework. 

Software versions: Requires at least [BEAST 2](https://github.com/CompEvol/beast2) v2.6.5 and [BEAST Labs](https://github.com/BEAST2-Dev/BEASTLabs) v1.9.0

The current release has the following features:

*Error Models*
* GT16 phased diploid nucleotide error model (16 states)
* Nucleotide error model (4 states - A, G, C, T)
* General error model (n states)
* Binary error model

*Substitution Models*
* GT16 phased diploid nucleotide model (16 states)
* Binary substitution model
* SiFit2 substitution model (2 states)
* SiFit3 substitution model (3 states)

## How to run

Download the executable jar file `phylonco.vx.x.x.jar`.

Execute the jar, and set the input file to `filepath/examples/*.xml`.

## How to build

Build BEAST 2 using `ant build_all_BEAST` for Linux/Mac or `ant windows` for Windows.

Then build BEAST Labs and Phylonco using `ant addon`.

## References
1. Bouckaert at al. (2019). [BEAST 2.5: An advanced software platform for Bayesian evolutionary analysis.](https://doi.org/10.1371/journal.pcbi.1006650)

2. Kozlov et al. (2021). [CellPhy: accurate and fast probabilistic inference of single-cell phylogenies from scDNA-seq data.](https://doi.org/10.1101/2020.07.31.230292)
 
3. Zafar et al. (2017). [SiFit: inferring tumor trees from single-cell sequencing data under finite-sites models.](https://doi.org/10.1186/s13059-017-1311-2)

4. Hou et al. (2012). [Single-cell exome sequencing and monoclonal evolution of a JAK2-negative myeloproliferative neoplasm.]( https://doi.org/10.1016/j.cell.2012.02.028)


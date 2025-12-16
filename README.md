# BEAST phylonco
This is a BEAST2 package for Bayesian inference of molecular data for cancer evolution. This package implements error models and substitution models for inference of timed trees in a Bayesian MCMC framework. 

[Paper](https://doi.org/10.1093/molbev/msac143): Chen, K., Moravec, J. C., Gavryushkin, A., Welch, D., & Drummond, A. J. (2022). Accounting for errors in data improves divergence time estimates in single-cell cancer evolution. Molecular biology and evolution, 39(8), msac143.

See [bioDS/beast-phylonco-paper](https://github.com/bioDS/beast-phylonco-paper) for datasets and analyses in the paper. 

If you have any questions, please use the `beast-users` google group https://groups.google.com/u/0/g/beast-users

## Software requirements

This package requires Java 17 and at least [BEAST v2.7](https://github.com/CompEvol/beast2)

For dependencies and lphy / lphybeast setup instructions see the [latest release notes](https://github.com/bioDS/beast-phylonco/releases/latest)

## Features
Phylonco versions v1.0.0 and above are able to use Beagle when calculating the tree likelihood with error models. 

**Beagle usage**: See instructions in [#running-beast-on-command-line](https://github.com/bioDS/beast-phylonco?tab=readme-ov-file#running-beast-on-command-line) below. 

* Using Beagle with GPU's can significantly improve computational speed by up to 10x compared to the standard Java version. This is important for larger datasets and models with larger matrices (e.g., GT16/GT10 models).  

The current release has the following features.

Error Models
* GT16 diploid nucleotide error model (16 states)
* Binary error model
* Nucleotide error model (4 states - A, G, C, T)
* General error model (n states)

Substitution Models
* GT16 diploid nucleotide model (16 states)
* Binary substitution model
* SiFit2 substitution model (2 states)
* SiFit3 substitution model (3 states)

## User guide
### How to install
You can install the Phylonco package using the Package Manager, or by command line.

**Package manager:**

1. Start `BEAUti`, open the `Package Manager` by selecting `File -> Manage packages` from the Menu.
2. From the `Package Manager`, select `phylonco` and click Install/Upgrade to install.
3. This may take a few minutes to install the package and dependencies. 
4. Restart `BEAUti` after installation is complete, and the `phylonco` package should now appear as Installed.

**Command line:**

```
~/beast/bin/packagemanager -add phylonco
```

### Input format
All models accept input genotypes in Nexus format. See [here](https://github.com/bioDS/beast-phylonco/blob/master/genotype_codes.pdf) for the genotype codes.

For the GT16 model, VCF files can be converted to Nexus or FASTA format using www.github.com/bioDS/vcf2fasta

### Running BEAST 
Start the BEAST software

Set the input file to one of the examples e.g., `examples/test_GT16_error.xml`

### Running BEAST on command line
These models can be run using the Java implementation (slower), or the Beagle implementation (faster). 
We recommend using [Beagle](https://github.com/beagle-dev/beagle-lib) with the GPU option for large datasets.

**Java option**

To run the Java implementation, use the commands:
```
Windows
java -jar c:\Users\BEASTUser\Desktop\BEAST\lib\launcher.jar -java beast.xml

Mac
/Applications/BEAST\ 2.6.6/bin/beast -java beast.xml

Linux
~/beast/bin/beast -java beast.xml
```

**Beagle option**

To run datasets using Beagle, we need to install [Beagle](https://github.com/beagle-dev/beagle-lib), and modify a single line in our XML. 
Open your XML in a text editor, find the element with `id='treeLikelihood'` and replace the spec field with `spec='phylonco.beast.evolution.likelihood.BeagleTreeLikelihoodWithError'`.

For example, replace this
```
<input id='treeLikelihood' spec='phylonco.beast.evolution.likelihood.TreeLikelihoodWithErrorFast' useAmbiguities='true' useTipLikelihoods='true'>
    ...
</input>
```

with the following
```
<input id='treeLikelihood' spec='phylonco.beast.evolution.likelihood.BeagleTreeLikelihoodWithError' useAmbiguities='true' useTipLikelihoods='true'>
    ...
</input>
```

To run Beagle with GPU use:
```
~/beast/bin/beast -beagle_GPU beast.xml
```

### How to setup model via BEAUti
BEAUti currently supports the GT16 Substitution model and Binary Substitution model. Support for error models in Beauti will be added soon. 

**Setting up a basic analysis**: 

1. Convert your data from VCF to Nexus format using the full diploid genotype option https://github.com/bioDS/vcf2fasta#converting-vcf-to-nexus
2. Open Beauti, from the menu go to `File -> Import Alignment` and select your Nexus file. Select the `Ambiguities` checkbox next to your alignment. 
3. Go to the `Site Model` panel, and select `GT16SubstitutionModel`.
4. To adjust tree model and priors, use the `Priors` panel.
5. Go to the `MCMC` panel and input your required chain length.
6. From the menu go to `File -> Save` to save your XML. This XML can now be run using BEAST 2.

**Adding error models to the analysis:**

This requires editing the XML generated from Beauti, see examples in the `examples` directory. 

## For developers 
Project now uses Maven, for details see [dev notes](https://github.com/LinguaPhylo/linguaPhylo/blob/master/DEV_NOTE.md)

For manual beast package installation instructions, see [here](http://www.beast2.org/managing-packages/)

## Citations
* BEAST v2.5: [Bouckaert at al. (2019)](https://doi.org/10.1371/journal.pcbi.1006650)

* BEAST2 Error models: [Chen et al. (2022)](https://doi.org/10.1093/molbev/msac143)

* GT16 model: [Kozlov et al. (2022)](https://doi.org/10.1186/s13059-021-02583-w) 
 
* SiFit2, SiFit3 models: [Zafar et al. (2017)](https://doi.org/10.1186/s13059-017-1311-2)

* Hou 2012 dataset: [Hou et al. (2012)]( https://doi.org/10.1016/j.cell.2012.02.028)


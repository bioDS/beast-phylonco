# BEAST Phylonco
This is a BEAST2 package for Bayesian inference of molecular data for cancer evolution. This package implements error models and substitution models for inference of timed trees in a Bayesian MCMC framework. 

[Paper](https://doi.org/10.1101/2021.03.17.435906): Chen et al. "Accounting for errors in data improves timing in single-cell cancer evolution." (2022)

See [bioDS/beast-phylonco-paper](https://github.com/bioDS/beast-phylonco-paper) for datasets and analyses in the paper. 

## Software requirements

This package requires Java 8, and at least
* [BEAST 2](https://github.com/CompEvol/beast2) v2.6.6 
* [BEAST Labs](https://github.com/BEAST2-Dev/BEASTLabs) v1.9.7

See dependencies in [version.xml](phylonco-lphybeast/version.xml)

## Features

The current release has the following features.

Error Models
* GT16 diploid nucleotide error model (16 states)
* Nucleotide error model (4 states - A, G, C, T)
* General error model (n states)
* Binary error model

Substitution Models
* GT16 diploid nucleotide model (16 states)
* Binary substitution model
* SiFit2 substitution model (2 states)
* SiFit3 substitution model (3 states)

## User guide
### How to install
You may install the Phylonco package using the Package Manager.

Start **BEAUti**, open the **Package Manager** by selecting `File -> Manage packages` from the Menu.

Click `Package repositories` to open a new popup window
<img src="https://raw.githubusercontent.com/rbouckaert/obama/master/doc/package_repos.png">

Click `Add URL` and add "https://raw.githubusercontent.com/CompEvol/CBAN/master/packages-extra.xml" to the entry.

Click `Done`. The Phylonco package should now appear in the **Package Manager**

From the **Package Manager**, select Phylonco and click Install/Upgrade to install.

### Input format
All models accept input genotypes in FASTA and Nexus formats. See [here](https://github.com/bioDS/beast-phylonco/blob/master/genotype_codes.pdf) for the genotype codes.

For the GT16 model, VCF files can be converted to FASTA format using www.github.com/bioDS/vcf2fasta

### Running BEAST 
Start the BEAST software

Make sure the "Use BEAGLE library" box is unchecked

Set the input file to one of the examples e.g., `examples/test_GT16_error.xml`


### Running BEAST on command line
Launch beast in java only mode by adding the `-java` option, e.g. 

```
Windows
java -jar c:\Users\BEASTUser\Desktop\BEAST\lib\launcher.jar -java beast.xml

Mac
/Applications/BEAST\ 2.6.6/bin/beast -java beast.xml

Linux
~/beast/bin/beast -java beast.xml
```

### How to run BEAUti
BEAUti UI is currently only supported for the Binary Substitution model. We are currently working on this :) 

Tutorials and new models will be added once Beauti support is completed.

## For developers 
To build use `./gradlew clean build` in the root directory, for details see [dev notes](https://github.com/LinguaPhylo/linguaPhylo/blob/master/DEV_NOTE.md)

For manual package installation instructions, see [here](http://www.beast2.org/managing-packages/)

## Citations
* BEAST v2.5: [Bouckaert at al. (2019)](https://doi.org/10.1371/journal.pcbi.1006650)

* BEAST2 Error models: [Chen et al. (2022)](https://doi.org/10.1101/2021.03.17.435906)

* GT16: [Kozlov et al. (2022)](https://doi.org/10.1186/s13059-021-02583-w) 
 
* SiFit2, SiFit3: [Zafar et al. (2017)](https://doi.org/10.1186/s13059-017-1311-2)

* Hou 2012 dataset: [Hou et al. (2012)]( https://doi.org/10.1016/j.cell.2012.02.028)


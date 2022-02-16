# BEAST Phylonco
This is a BEAST2 package for Bayesian inference of molecular data for cancer evolution. This package implements error models and substitution models for inference of timed trees in a Bayesian MCMC framework. 

## Software requirements

This package requires Java 8, and at least
* [BEAST 2](https://github.com/CompEvol/beast2) v2.6.5 
* [BEAST Labs](https://github.com/BEAST2-Dev/BEASTLabs) v1.9.0

## Features

The current release has the following features.

Error Models
* GT16 phased diploid nucleotide error model (16 states)
* Nucleotide error model (4 states - A, G, C, T)
* General error model (n states)
* Binary error model

Substitution Models
* GT16 phased diploid nucleotide model (16 states)
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

### Running BEAST 
Start the BEAST software

Make sure the "use BEAGLE" box is unchecked

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

### How to run Beauti
Beauti UI is currently not supported. We are currently working on this :) 

Tutorials and new models will be added once Beauti support is completed.

## For developers 
### How to build 

Build BEAST 2 using `ant build_all_BEAST` for Linux/Mac or `ant windows` for Windows.

Then build BEAST Labs and Phylonco using `ant addon`.

### Manual install 
See [here](http://www.beast2.org/managing-packages/) for manual package install instructions.

## References
1. Bouckaert at al. (2019). [BEAST 2.5: An advanced software platform for Bayesian evolutionary analysis.](https://doi.org/10.1371/journal.pcbi.1006650)

2. Kozlov et al. (2021). [CellPhy: accurate and fast probabilistic inference of single-cell phylogenies from scDNA-seq data.](https://doi.org/10.1101/2020.07.31.230292)
 
3. Zafar et al. (2017). [SiFit: inferring tumor trees from single-cell sequencing data under finite-sites models.](https://doi.org/10.1186/s13059-017-1311-2)

4. Hou et al. (2012). [Single-cell exome sequencing and monoclonal evolution of a JAK2-negative myeloproliferative neoplasm.]( https://doi.org/10.1016/j.cell.2012.02.028)


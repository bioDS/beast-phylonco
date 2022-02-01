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

*For scripting support, see [bioDS/lphy-phylonco](https://github.com/bioDS/lphy-phylonco)

## How to install
You may install the Phylonco package either using the Package Manager, or manually.

### Package manager install
Start **BEAUti**, open the **Package Manager** by selecting `File -> Manage packages` from the Menu.

Click `Package repositories` to open a new popup window
<img src="https://raw.githubusercontent.com/rbouckaert/obama/master/doc/package_repos.png">

Click `Add URL` and add "https://raw.githubusercontent.com/CompEvol/CBAN/master/packages-extra.xml" to the entry.

Click `Done`. The Phylonco package should now appear in the **Package Manager**

From the **Package Manager**, select Phylonco and click Install/Upgrade to install.

---

### Manual install
Download the latest addon zip `phylonco.addon.vx.x.x.zip` from the [releases page](https://github.com/bioDS/beast-phylonco/releases/)

Create a new subdirectory `Phylonco` at the location: 
```
For Windows in Users\<username>\BEAST\2.x\Phylonco
For Mac in /Users/<username>/Library/Application Support/BEAST/2.x/Phylonco
For Linux in /home/<username>/.beast/2.x/Phylonco
```
Here `<username>` is your username, “2.x” refers to the major version of BEAST, for example 2.x = 2.1 for BEAST version 2.1.3.
 
Then extract the contents of the addon zip to the `Phylonco` subdirectory.

For BEAST v2.5.x and later, you need to reset the class path stored in the beauti.properties file. The easiest way to do this is by starting BEAUti and selecting `File -> Clear class path` from the Menu.

See [here](http://www.beast2.org/managing-packages/) for general package install instructions.

## How to run 

Start **BEAST2**, and set the input file to one of the examples in `examples/*.xml`.

## User interface 
Beauti UI is currently not supported. We are currently working on this :) 

Tutorials and new models will be added once Beauti support is completed.

## How to build (developers)

Build BEAST 2 using `ant build_all_BEAST` for Linux/Mac or `ant windows` for Windows.

Then build BEAST Labs and Phylonco using `ant addon`.

## References
1. Bouckaert at al. (2019). [BEAST 2.5: An advanced software platform for Bayesian evolutionary analysis.](https://doi.org/10.1371/journal.pcbi.1006650)

2. Kozlov et al. (2021). [CellPhy: accurate and fast probabilistic inference of single-cell phylogenies from scDNA-seq data.](https://doi.org/10.1101/2020.07.31.230292)
 
3. Zafar et al. (2017). [SiFit: inferring tumor trees from single-cell sequencing data under finite-sites models.](https://doi.org/10.1186/s13059-017-1311-2)

4. Hou et al. (2012). [Single-cell exome sequencing and monoclonal evolution of a JAK2-negative myeloproliferative neoplasm.]( https://doi.org/10.1016/j.cell.2012.02.028)


#!/bin/bash
# Description: sequence simulation based on SiFit model by Zafar et al. (2017)
# Usage: ./run_simulation.sh
mkdir -p output
mkdir -p output/sequences
mkdir -p output/xml
cd r

echo "installing r packages"
Rscript -e 'install.packages("expm", repos="https://cran.rstudio.com")'

echo "running simulation"
Rscript "simulate_tree.r"

cd ..
echo "done!"

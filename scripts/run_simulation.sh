#!/bin/bash
# Description: sequence simulation based on SiFit model by Zafar et al. (2017)
# Usage: ./run_simulation.sh

echo "installing r packages"
Rscript -e 'install.packages("expm", repos="https://cran.rstudio.com")'

echo "running simulation"
cd r
Rscript "simulate_tree.r"

cd ..
echo "done!"

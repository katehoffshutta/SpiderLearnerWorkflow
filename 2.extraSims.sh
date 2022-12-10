#!/bin/bash

# This file runs all the additional simulations outside of the four main settings (ABCD)
# in alphabetical order.
# Note that each R script relies on parameters in config.yml. In particular, 
# the configs are set to run a "pilot" analysis that has
# trivial amount of simulations (2), folds (2), and cores (1-2). 
# Adjust as needed to reproduce manuscript results and to match computing power.

# bootstrap_comparison
Rscript --vanilla Simulations/comparisonWithBootstrap.R > Logs/bootstrap_out.txt

# bounded_loss
Rscript --vanilla Simulations/boundedLossStudy.R > Logs/bounded_loss_out.txt

# kfold (choice of k)
Rscript --vanilla Simulations/kfoldStudy.R > Logs/kfold_out.txt

# library_sensitivity
Rscript --vanilla Simulations/librarySensitivityStudy.R > Logs/library_out.txt

# mle_study
Rscript --vanilla Simulations/mleStudy.R > Logs/mle_out.txt

# numerical_stability
Rscript --vanilla Simulations/numericalOptimizationStudy.R > Logs/numerical_stability_out.txt

# ovarian_library_study
Rscript --vanilla Simulations/ovarianLibraryStudy.R > Logs/ovarian_out.txt

# runtime_study
Rscript --vanilla Simulations/runtimeStudy.R > Logs/runtime_out.txt



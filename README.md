# README.md

# SpiderLearner Workflow

This repository contains the code for simulations testing the performance of our R package `ensembleGGM` and an illustrative application to ovarian cancer data. 

Ensemble learning can be computationally expensive, and so most of the code in this repository is designed to be run in a large-scale computing environment with multiple cores. Our simulations were conducted primarily on the Massachusetts Green High Performance Computing Cluster (https://www.mghpcc.org/). 

To allow the user to pilot these simulations and applications locally and allocate resources depending on what they have available, we have used the `config` R package along with configuration files `Simulations/config.yml` and `Applications/config.yml` to encode pilot versions of the scripts that include a small number of folds (2-3) and a small number of cores (1-2). We recommend using K=10 folds in general. The number of cores can be adjusted as desired based on the available hardware.  

Please use Github issues to report general questions. For individual-specific questions, email Kate Hoff Shutta at kshutta at hsph.harvard.edu.



#!/bin/bash
#BSUB -n 40
#BSUB -W 27:59
#BSUB -R rusage[mem=2000]
#BSUB -R span[hosts=1]
module load R/3.6.0
module load gcc/8.1.0
Rscript --vanilla simulations_C.R > c.txt

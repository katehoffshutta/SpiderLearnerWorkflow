#!/bin/bash
#BSUB -n 10
#BSUB -W 03:59
#BSUB -R rusage[mem=3000]
#BSUB -R span[hosts=1]
module load R/3.6.0
module load gcc/8.1.0
Rscript --vanilla ovarianSmall.R > ovarian_out.txt

# packages from CRAN
install.packages(clime)
install.packages(config)
install.packages(devtools)
install.packages(doParallel)
install.packages(dplyr)
install.packages(ensembleGGM)
install.packages(foreach)
install.packages(ggplot2)
install.packages(huge)
install.packages(igraph)
install.packages(MASS)
install.packages(MatrixCorrelation)
install.packages(moments)
install.packages(pracma)

# BiocManager if needed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# packages from Bioconductor
BiocManager::install("affy")
BiocManager::install("curatedOvarianData")

# ensembleGGM from github
install_github("katehoffshutta/ensembleGGM")
#  This script generates the file mxDist.csv and plots the partial correlation distribution for Supplementary Figure 1 and plots the structures for Figure 3
#  The raw data are not provided; the purpose of this script is to show the reader how we made the distribution.
# The file mxDist.csv is provided in the repository.

library(ggplot2)
library(huge)
library(igraph)
library(MASS)
library(moments)
library(pracma)

#### Get discrete uniform distribution from CATHGEN Example ####
standardize = function(x) return((x-mean(x))/sd(x))

cathgen = as.matrix(read.table("data.txt",sep="\t", header=T))
cathgen=apply(cathgen,2,standardize)
cathgenGGM = huge(cathgen, method="glasso")
cathgenOptGGM = huge.select(cathgenGGM, criterion="ebic",ebic.gamma=0)
cathgenEdgeWeights = c(-cov2cor(cathgenOptGGM$opt.icov)+diag(ncol(cathgen)))
cathgenNonZeroEdgeWeights = cathgenEdgeWeights[cathgenEdgeWeights != 0]
cathgenHist = hist(cathgenNonZeroEdgeWeights,breaks=100)
write.csv(data.frame("mids"=cathgenHist$mids, "density"=cathgenHist$density[1:108]),file="mxDist.csv",row.names=F,quote=F)

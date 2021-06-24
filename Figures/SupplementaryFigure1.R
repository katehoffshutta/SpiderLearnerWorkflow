library(huge)
library(igraph)
library(moments)
library(pracma)

#### Get discrete uniform distribution from CATHGEN Example ####
standardize = function(x) return((x-mean(x))/sd(x))

cathgen = as.matrix(read.table("~/umass/research/networks/metabolomics2019/Biostat-study-master/Data/data.txt",sep="\t", header=T))
cathgen=apply(cathgen,2,standardize) # gave same results - glasso standardizes
cathgenGGM = huge(cathgen, method="glasso")
cathgenOptGGM = huge.select(cathgenGGM, criterion="ebic",ebic.gamma=0)
cathgenEdgeWeights = c(-cov2cor(cathgenOptGGM$opt.icov)+diag(ncol(cathgen)))
cathgenNonZeroEdgeWeights = cathgenEdgeWeights[cathgenEdgeWeights != 0]
cathgenHist = hist(cathgenNonZeroEdgeWeights,breaks=100)

pdf("cathgenPrecisionDist.pdf",width=8,height=8)
hist((cathgenOptGGM$opt.icov- diag(cathgenOptGGM$opt.icov))[which(cathgenOptGGM$opt.icov- diag(cathgenOptGGM$opt.icov) != 0)],
     breaks=100,xlab="Precision Matrix Entries", main="CATHGEN Off-diagonal Precision Matrix Distribution")
dev.off()

skewness(cathgenNonZeroEdgeWeights)
kurtosis(cathgenNonZeroEdgeWeights)


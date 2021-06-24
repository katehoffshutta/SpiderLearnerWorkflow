library(huge)
library(igraph)
library(MASS)
library(pracma)

source("errorMetrics.R")
source("spiderLearnerOOP.R")
source("spiderLearnerCandidates.R")

standardize = function(x){return((x-mean(x))/sd(x))}
sampleNetworkData = function(N, covMat)
{
  sample = mvrnorm(N, mu = rep(0,nrow(covMat)), covMat)
  sample = apply(sample,2,standardize)
  return(sample)
}

boost = function(myMatrix)
{
  minEig = min(eigen(myMatrix)$values)
  if(minEig < 0)
  {
    print("Boosting Needed!")
    return(myMatrix - minEig*1.01*diag(ncol(myMatrix)))
  }
  else
    return(myMatrix)
}

cathgenHist = read.csv("mxDist.csv")

#### Make Gold Standard Network Structures #### 

set.seed(42)
makeGoldStandardNetsB = function(P)
{
  PminusOne = P-1
  
  # make random graph
  er1 = sample_gnp(P,0.05)
  er2 = sample_gnp(P,0.1)
  er3 = sample_gnp(P,0.25)
  er4 = sample_gnp(P,0.5)
  er5 = sample_gnp(P,0.75)
  er6 = sample_gnp(P,1)

  
  #### Assign Edge Weights ####
  
  assignRunifEdgeWeights = function(graph, edgeMin=-1, edgeMax=1)
  {
    P = length(V(graph))
    weights = sample(c(-1,1),size=length(E(graph)),replace=T)*runif(n = length(E(graph)), min = edgeMin, max = edgeMax) 
    E(graph)$weights = weights
    precMat = -1*as_adjacency_matrix(graph, attr = c("weights"), type="both") + diag(P)
    return(list(graph, precMat))
  }
  
  assignCathgenEdgeWeights = function(graph,ehist)
  {
    P = length(V(graph))
    weights = sample(ehist$mids, replace=T,size=length(E(graph)), prob=ehist$density)  
    E(graph)$weights = weights
    precMat = -1*as_adjacency_matrix(graph, attr = c("weights"), type="both") + diag(P)
    return(list(graph, precMat))
  }
  
  
  er1Cathgen = assignCathgenEdgeWeights(er1,cathgenHist)
  er2Cathgen = assignCathgenEdgeWeights(er2,cathgenHist)
  er3Cathgen = assignCathgenEdgeWeights(er3,cathgenHist)
  er4Cathgen = assignCathgenEdgeWeights(er4,cathgenHist)
  er5Cathgen = assignCathgenEdgeWeights(er5,cathgenHist)
  er6Cathgen = assignCathgenEdgeWeights(er6,cathgenHist)

  adjMatListCathgen = c(boost(er1Cathgen[[2]]),
                        boost(er2Cathgen[[2]]),
                        boost(er3Cathgen[[2]]),
                        boost(er4Cathgen[[2]]),
                        boost(er5Cathgen[[2]]),
                        boost(er6Cathgen[[2]]))
  for(mat in adjMatListCathgen)
  {
    print(isposdef(as.matrix(mat)))
    print(min(eigen(mat)$values))
  }
  
  return(adjMatListCathgen)
}

nPred = 50
nObs = 1600
set.seed(42)

sparseAdjMats = makeGoldStandardNetsB(nPred)

s = SpiderLearner$new()

apple = HugeEBICCandidate$new(gamma = 0)
banana = HugeEBICCandidate$new(gamma = 0.5)
clementine = HugeRICCandidate$new()
date = HGlassoCandidate$new()
elderberry = MLECandidate$new()
fraise = HugeStARSCandidate$new(thres = 0.05)
grape = HugeStARSCandidate$new(thres = 0.1)
honeydew = QGraphEBICCandidate$new(gamma = 0)
icewine = QGraphEBICCandidate$new(gamma = 0.5)

candidates = list(apple, 
                  banana, 
                  clementine, 
                  date, 
                  elderberry,
                  fraise,
                  grape,
                  honeydew,
                  icewine)

for(candidate in candidates)
{
  s$addCandidate(candidate)
}

nIt = 30
ensModelWeights = array(rep(NA,nIt*length(sparseAdjMats)*length(candidates)),dim=c(nIt,length(sparseAdjMats),length(candidates)))

for(N in 1:nIt)
{
  print(paste("Iteration",N))
  for(i in 1:length(sparseAdjMats))
  {
    print(paste("Sparse adj mat",i))
    thisNetwork = sparseAdjMats[[i]]
    thisSample = sampleNetworkData(N=nObs,covMat=solve(thisNetwork))
    thisTestSample = sampleNetworkData(N=nObs,covMat=solve(thisNetwork))
    ensModel = s$runSpiderLearner(thisSample,10,standardize=FALSE,nCores = 10)
    ensModelWeights[N,i,] = ensModel$weights
  }
}

save(ensModelWeights,file="mleWeights.rda")

# 
# sparsityVec = c(0.05,0.10,0.25,0.5,0.75,1)
# pdf("mleWeight_bigUnif.pdf",width=8,height=8)
# plot(sparsityVec,ensModelWeights[,8],ylim=c(0,1),xlab="density = 1-sparsity",ylab="weight in ensemble model",pch=20,
#      main="Weight of MLE in Ensemble \n Depends on Sparsity")
# lines(sparsityVec,ensModelWeights[,8],lty=1,lwd=2)
# for(i in 1:7)
# {
#   points(sparsityVec,ensModelWeights[,i],lty=2,col=rainbow(7)[i],pch=20)
#   lines(sparsityVec,ensModelWeights[,i],lty=2,col=rainbow(7)[i],lwd=1)
# }
# legend(0.6,0.6,c(names(erHigh_RFN)[1:7],"MLE"),col=c(rainbow(7),"black"),pch=20,lty=c(rep(2,7),1),cex=0.8)
# dev.off()

library(config)
library(ggplot2)
library(huge)
library(igraph)
library(MASS)
library(moments)
library(pracma)
library(ensembleGGM)

config = config::get(config="pilot_sim_bounded_loss", file="Simulations/config.yml")

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

#### Get discrete uniform distribution from REALDATA Example ####

realDataHist = read.csv("Simulations/mxDist.csv")

#### Make Gold Standard Network Structures #### 

makeGoldStandardNet = function(P)
{
  PminusOne = P-1
  
  # make random graph,high Density
  er1 = sample_gnp(P,0.20)

  
  #### Assign Edge Weights ####
 
  
  assignRealDataEdgeWeights = function(graph,ehist)
  {
    P = length(V(graph))
    weights = sample(ehist$mids, replace=T,size=length(E(graph)), prob=ehist$density)  
    E(graph)$weights = weights
    precMat = -1*as_adjacency_matrix(graph, attr = c("weights"), type="both") + diag(P)
    return(list(graph, precMat))
  }
  
  
  er1RealData = assignRealDataEdgeWeights(er1,realDataHist)
 
  adjMatRealData = boost(er1RealData[[2]])

  print(isposdef(as.matrix(adjMatRealData)))
 
  return(adjMatRealData)
}

nPred = 50
nObs = 1600
set.seed(42)

adjMat = makeGoldStandardNet(nPred)

s = SpiderLearner$new()

apple = HugeEBICCandidate$new(gamma = 0)
banana = HugeEBICCandidate$new(gamma = 0.5)
clementine = HugeRICCandidate$new()
date = HGlassoCandidate$new()
elderberry = MLECandidate$new()

candidates = list(apple, 
                  banana, 
                  clementine, 
                  date, 
                  elderberry)

for(candidate in candidates)
{
  s$addCandidate(candidate)
}

nIt = config$nSim
nFolds = config$nFolds

ensModelWeightsUnbd = matrix(rep(NA,nIt*length(candidates)),ncol=length(candidates))
ensModelWeightsBd = matrix(rep(NA,nIt*length(candidates)),ncol=length(candidates))

adjCov = solve(adjMat)

for(N in 1:nIt)
{
  print(paste("Iteration",N))

  thisSample = sampleNetworkData(N=nObs,covMat=adjCov)
  thisTestSample = sampleNetworkData(N=nObs,covMat=adjCov)
  ensModel = s$runSpiderLearner(thisSample,nFolds,standardize=FALSE,boundedLoss= F)
  ensModelWeightsUnbd[N,] = ensModel$weights
  ensModel = s$runSpiderLearner(thisSample,nFolds,standardize=FALSE,boundedLoss= T)
  ensModelWeightsBd[N,] = ensModel$weights
}

write.csv(ensModelWeightsUnbd,file="Results/Pilot/ensModelWeightsUnbd.csv",row.names=F,quote=F)
write.csv(ensModelWeightsBd,file="Results/Pilot/ensModelWeightsBd.csv",row.names=F,quote=F)

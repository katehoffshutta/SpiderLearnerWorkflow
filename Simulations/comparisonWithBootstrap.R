library(config)
library(dplyr)
library(ensembleGGM)
library(MASS)
library(huge)

source("Simulations/errorMetrics.R")

args = commandArgs(trailingOnly=TRUE)
config_name = args[1]

# read in config parameters from config file
config = config::get(config="pilot_sim_bootstrap", file="Simulations/config.yml")

load("Results/eightNetworksAC.rda")
eightNetworks = ac
precMat = eightNetworks[[2]] # random graph high density
print(nrow(precMat))
# plot(graph_from_adjacency_matrix(-cov2cor(as.matrix(precMat)),weighted=T,diag=F,mode="undirected"),layout=layout_with_lgl)

standardize = function(x){return((x-mean(x))/sd(x))}
sampleNetworkData = function(N, covMat)
{
  sample = mvrnorm(N, mu = rep(0,nrow(covMat)), covMat)
  sample = apply(sample,2,standardize)
  return(sample)
}


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

s = MakeSpiderLearner(candidates)

nObs = 1000
nIt = config$nSim

compResults = list()

set.seed(42)
for(i in 1:nIt)
{ 
  print(paste("[comparisonWithBootstrap] Iteration:", i))
  thisSample = sampleNetworkData(N=nObs,covMat=solve(precMat))
  
  ensModel = s$runSpiderLearner(thisSample,config$nFolds,
                                standardize=FALSE,
                                boundedLoss= F,
                                nCores = config$nCores)
  
  # now test bootstrap approach
  bootstrapModels = list()
  B = config$nBoot

  for(b in 1:B)
  {
    bootIdx = sample(1:nrow(thisSample),replace=T)
    bootSample = thisSample[bootIdx,]
    myMod = huge(bootSample,method = "glasso")
    myModOpt = huge.select(myMod, criterion = "ebic", ebic.gamma = 0)
    bootstrapModels[[b]] = myModOpt$opt.icov
  }
  
  bootEnsemble = 1/B*Reduce('+',bootstrapModels)
  
  # compare likelihood on new data
  testSample = sampleNetworkData(N=nObs,covMat=solve(precMat))
  bootstrapOOSL = loglikLossfunction(bootEnsemble,testSample)
  ensOOSL = loglikLossfunction(ensModel$optTheta,testSample)
  diffOOSL = ensOOSL - bootstrapOOSL # positive: ensemble outperforms bootstrap
  thisRes = data.frame("bootstrapOOSL"=bootstrapOOSL,
                 "ensOOSL"=ensOOSL,
                 "diffOOSL"= diffOOSL)
  
  compResults[[i]] = thisRes
  dfOut = bind_rows(compResults)
  write.csv(dfOut, "Results/Pilot/bootstrapComparisonResults.csv",quote=F,row.names=F)

}

